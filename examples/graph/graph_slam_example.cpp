//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//
// 3D pose-graph SLAM with odometry and pinhole-camera observations of 3D
// landmarks, via the GraphProblem interface.
//
// World model
// -----------
//   * N robot poses on SE(3) parametrised as (p, q) with
//        p = (px, py, pz) in R^3,
//        q = (qw, qx, qy, qz) a unit quaternion.
//     The block primal dimension is 7. The tangent dimension is 6:
//        delta = (delta_p ; delta_phi),  delta_phi in so(3) (right-trivialised).
//     Retraction:
//        p_next = p + alpha * delta_p,                  (world-frame translation)
//        q_next = q (x) Exp(alpha * delta_phi),         (right-multiplied so(3) exp)
//
//   * M static 3D landmarks l_j = (lx, ly, lz) in R^3. Primal = tangent = 3.
//
// Observations
// ------------
//   * Odometry: for each consecutive pair (k, k+1), a noisy measurement of
//     the body-frame relative pose:
//        delta_xy = R(q_k)^T (p_{k+1} - p_k),
//        delta_phi = log( R(q_k)^T R(q_{k+1}) ) in so(3).
//
//   * Pinhole landmark observations: at every pose k that "sees" landmark j
//     (i.e. it falls in front of the camera within VIEW_RANGE and inside the
//     image bounds), the camera emits a 2D pixel measurement (u, v) where
//        l_cam = R(q_k)^T (l_j - p_k),
//        u = f * l_cam_x / l_cam_z,
//        v = f * l_cam_y / l_cam_z.
//     Per-pose the K closest in-FOV landmarks are kept (MAX_OBS_PER_POSE) to
//     mimic the bounded feature-track count of real SLAM front-ends.
//
// Scenario
// --------
//   The synthetic data layer (slam_synthetic_data.hpp) builds an aerial
//   drone "lawn-mower" survey: the drone flies n_rows parallel rows at a
//   constant altitude with a fixed downward-looking gimballed camera and
//   observes landmarks scattered on the ground. Loop closures arise from
//   landmarks that fall in the small sidelap zone between adjacent rows.
//
// Code structure
// --------------
//   1. SyntheticSlamData / generate_synthetic_data: pure scenario generation —
//      ground truth, noisy measurements, dead-reckoning initial guess. Knows
//      nothing about fatrop, GraphAbstract or NLP variables.
//   2. SlamProblem: takes a SlamProblemInput (prior + measurements + initial
//      guess; no ground truth) and implements GraphAbstract — it exposes the
//      block sizes, the Hessian / gradient / Jacobian callbacks and the SE(3)
//      retraction.
//

#include <fatrop/fatrop.hpp>
#include <fatrop/graph/graph_abstract.hpp>
#include <fatrop/graph/hessian.hpp>
#include <fatrop/graph/jacobian.hpp>
#include <fatrop/graph/nlp_graph.hpp>
#include <fatrop/graph/problem_info.hpp>
#include <fatrop/graph/problem_type.hpp>
#include <fatrop/ip_algorithm/ip_alg_builder.hpp>
#include <fatrop/ip_algorithm/ip_algorithm.hpp>
#include <fatrop/ip_algorithm/ip_data.hpp>

#include "slam_helpers.hpp"
#include "slam_synthetic_data.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <utility>
#include <vector>

using namespace fatrop;

// ============================================================================
// LAYER 2 — Problem formulation (GraphAbstract)
//   - Pure NLP / fatrop wiring: block sizes, retraction, cost, gradient,
//     Hessian, inequality constraints. Knows nothing about ground truth or
//     about the synthetic data generator — driven entirely by the neutral
//     SlamProblemInput below.
// ============================================================================

// Bring the SE(3) / camera types and helpers (quat_to_rot, project,
// rel_pose_jacobian, ...) plus OdomMeas / PixelMeas into scope. See
// slam_helpers.hpp.
using namespace fatrop::examples;

// Everything the SlamProblem needs to formulate the bundle-adjustment NLP.
// The driver populates this from whatever source it likes (synthetic data,
// a recorded log, a live front-end).
struct SlamProblemInput
{
    Pose prior;                          // anchor prior on pose 0
    std::vector<OdomMeas> odom;          // body-frame relative-pose measurements
    std::vector<PixelMeas> observations; // pinhole landmark observations
    std::vector<Pose> pose_init;         // initial guess for the poses (size = n_pose)
    std::vector<Landmark> lmk_init;      // initial guess for the landmarks (size = n_landmark)
};

struct SlamProblemParams
{
    Scalar prior_info = 1e4;
    Scalar odom_info = 1e2;
    Scalar pix_info = 1.0;
};

class SlamProblem : public GraphAbstract
{
public:
    SlamProblem(SlamProblemInput input, const SlamProblemParams &p = {})
        : in_(std::move(input)), p_(p),
          n_pose_(static_cast<Index>(in_.pose_init.size())),
          n_landmark_(static_cast<Index>(in_.lmk_init.size()))
    {
        build_edge_indices();
        // Scale the objective by the total number of scalar residuals so
        // its magnitude stays O(1) regardless of the dataset size. Without
        // this the initial cost is on the order of 1e6+, which makes the
        // IP's natural tolerances (absolute, scale-invariant) feel
        // unbalanced relative to the per-residual contributions.
        const Index n_res = 6 + 6 * static_cast<Index>(in_.odom.size()) +
                            2 * static_cast<Index>(in_.observations.size());
        cost_scale_ = 1.0 / static_cast<Scalar>(n_res);
    }

    Index n_pose() const { return n_pose_; }
    Index n_landmark() const { return n_landmark_; }

    // Primal-vector layout: [pose_0 | ... | pose_{n_pose-1} | lmk_0 | ... | lmk_{n_lmk-1}].
    Index pose_primal_offset(Index k) const { return k * POSE_PRIM_DIM; }
    Index lmk_primal_offset(Index j) const { return n_pose_ * POSE_PRIM_DIM + j * LMK_DIM; }

    // ---- Block dimensions --------------------------------------------------
    Index get_num_blocks() const override { return n_pose_ + n_landmark_; }
    Index get_block_size(Index k) const override
    {
        return (k < n_pose_) ? POSE_PRIM_DIM : LMK_DIM;
    }
    Index get_block_size_tan(Index k) const override
    {
        return (k < n_pose_) ? POSE_TAN_DIM : LMK_DIM;
    }

    std::vector<std::pair<Index, Index>> get_off_diag_edges() const override
    {
        std::vector<std::pair<Index, Index>> edges;
        for (Index k = 0; k + 1 < n_pose_; ++k)
            edges.emplace_back(pose_block(k), pose_block(k + 1));
        for (const auto &obs : in_.observations)
            edges.emplace_back(pose_block(obs.pose_k), lmk_block(obs.lmk_j));
        return edges;
    }

    // ---- Initial guess -----------------------------------------------------
    Index get_initial(Scalar *x) const override
    {
        for (Index k = 0; k < n_pose_; ++k)
            for (Index a = 0; a < POSE_PRIM_DIM; ++a)
                x[pose_primal_offset(k) + a] = in_.pose_init[k][a];
        for (Index j = 0; j < n_landmark_; ++j)
            for (Index a = 0; a < LMK_DIM; ++a)
                x[lmk_primal_offset(j) + a] = in_.lmk_init[j][a];
        return 0;
    }

    // ---- SE(3) per-block retraction ---------------------------------------
    void apply_retraction(const Scalar *x, const Scalar *delta_x, const Scalar alpha,
                          Scalar *x_next) override
    {
        for (Index k = 0; k < n_pose_; ++k)
        {
            const Index po = pose_primal_offset(k);
            const Index to = k * POSE_TAN_DIM;
            x_next[po + 0] = x[po + 0] + alpha * delta_x[to + 0];
            x_next[po + 1] = x[po + 1] + alpha * delta_x[to + 1];
            x_next[po + 2] = x[po + 2] + alpha * delta_x[to + 2];
            const Scalar phi[3] = {alpha * delta_x[to + 3], alpha * delta_x[to + 4],
                                   alpha * delta_x[to + 5]};
            Scalar dq[4];
            quat_exp(phi, dq);
            Scalar q_new[4];
            quat_mul(x + po + 3, dq, q_new);
            quat_normalize(q_new);
            for (Index a = 0; a < 4; ++a)
                x_next[po + 3 + a] = q_new[a];
        }
        for (Index j = 0; j < n_landmark_; ++j)
        {
            const Index lo = lmk_primal_offset(j);
            const Index to = n_pose_ * POSE_TAN_DIM + j * LMK_DIM;
            for (Index a = 0; a < LMK_DIM; ++a)
                x_next[lo + a] = x[lo + a] + alpha * delta_x[to + a];
        }
    }

    // ---- Objective ---------------------------------------------------------
    Index eval_f(const Scalar *objective_scale, const Scalar *x, Scalar *res) override
    {
        const Scalar s = objective_scale[0];
        Scalar v = 0.0;
        v += 0.5 * prior_w() * prior_sq_residual(x);
        for (const auto &meas : in_.odom)
        {
            Scalar r[6];
            odom_residual(x, meas, r);
            for (Index a = 0; a < 6; ++a)
                v += 0.5 * odom_w() * r[a] * r[a];
        }
        for (const auto &obs : in_.observations)
        {
            Pose pose;
            Landmark lmk;
            load_pose(x, obs.pose_k, pose);
            load_lmk(x, obs.lmk_j, lmk);
            Scalar pix[2];
            if (!project(pose, lmk, pix))
                continue;
            const Scalar ru = pix[0] - obs.uv[0];
            const Scalar rv = pix[1] - obs.uv[1];
            v += 0.5 * pix_w() * (ru * ru + rv * rv);
        }
        *res = s * v;
        return 0;
    }

    // ---- Gradient ----------------------------------------------------------
    Index eval_grad_k(Index k, const Scalar *objective_scale, const Scalar *x,
                      Scalar *res) override
    {
        const Scalar s = objective_scale[0];
        const Index n_tan = get_block_size_tan(k);
        for (Index a = 0; a < n_tan; ++a)
            res[a] = 0.0;
        if (k == 0)
        {
            Scalar gp[6];
            prior_gradient(x, gp);
            for (Index a = 0; a < 6; ++a)
                res[a] += s * prior_w() * gp[a];
        }
        if (k < n_pose_)
        {
            for (Index e_idx : pose_odom_edges_[k])
            {
                const OdomMeas &meas = in_.odom[e_idx];
                Scalar r[6];
                odom_residual(x, meas, r);
                Scalar Ja[6][6], Jb[6][6];
                odom_jacobian(x, meas, Ja, Jb);
                if (meas.a == k)
                    add_JT_r_pose(res, Ja, r, s * odom_w());
                if (meas.b == k)
                    add_JT_r_pose(res, Jb, r, s * odom_w());
            }
            for (Index e_idx : pose_obs_edges_[k])
            {
                const PixelMeas &obs = in_.observations[e_idx];
                Pose pose;
                Landmark lmk;
                load_pose(x, obs.pose_k, pose);
                load_lmk(x, obs.lmk_j, lmk);
                Scalar pix[2];
                if (!project(pose, lmk, pix))
                    continue;
                const Scalar r[2] = {pix[0] - obs.uv[0], pix[1] - obs.uv[1]};
                Scalar Jp[2][6], Jl[2][3];
                project_jacobian(pose, lmk, Jp, Jl);
                (void)Jl;
                for (Index a = 0; a < 6; ++a)
                    res[a] += s * pix_w() * (Jp[0][a] * r[0] + Jp[1][a] * r[1]);
            }
        }
        else
        {
            const Index j = k - n_pose_;
            for (Index e_idx : lmk_obs_edges_[j])
            {
                const PixelMeas &obs = in_.observations[e_idx];
                Pose pose;
                Landmark lmk;
                load_pose(x, obs.pose_k, pose);
                load_lmk(x, obs.lmk_j, lmk);
                Scalar pix[2];
                if (!project(pose, lmk, pix))
                    continue;
                const Scalar r[2] = {pix[0] - obs.uv[0], pix[1] - obs.uv[1]};
                Scalar Jp[2][6], Jl[2][3];
                project_jacobian(pose, lmk, Jp, Jl);
                (void)Jp;
                for (Index a = 0; a < 3; ++a)
                    res[a] += s * pix_w() * (Jl[0][a] * r[0] + Jl[1][a] * r[1]);
            }
        }
        return 0;
    }

    // ---- Hessian (Gauss-Newton) -------------------------------------------
    Index eval_Hk(Index i, Index j, const Scalar *objective_scale, const Scalar *x,
                  const Scalar * /*mult*/, MAT *res) override
    {
        const Scalar s = objective_scale[0];
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        if (i == 0 && j == 0)
            for (Index a = 0; a < 6; ++a)
                blasfeo_matel_wrap(res, a, a) += s * prior_w();

        const bool is_pp = (i < n_pose_ && j < n_pose_);
        const bool is_pl = (i >= n_pose_ && j < n_pose_);
        const bool is_lp = (i < n_pose_ && j >= n_pose_);
        const bool is_ll = (i >= n_pose_ && j >= n_pose_);

        if (is_pp)
        {
            if (i == j)
            {
                for (Index e_idx : pose_odom_edges_[i])
                {
                    const OdomMeas &meas = in_.odom[e_idx];
                    Scalar Ja[6][6], Jb[6][6];
                    odom_jacobian(x, meas, Ja, Jb);
                    if (meas.a == i)
                        add_JTJ_pose(res, Ja, Ja, s * odom_w());
                    if (meas.b == i)
                        add_JTJ_pose(res, Jb, Jb, s * odom_w());
                }
                for (Index e_idx : pose_obs_edges_[i])
                {
                    const PixelMeas &obs = in_.observations[e_idx];
                    Pose pose;
                    Landmark lmk;
                    load_pose(x, obs.pose_k, pose);
                    load_lmk(x, obs.lmk_j, lmk);
                    Scalar Jp[2][6], Jl[2][3];
                    if (!project_jacobian(pose, lmk, Jp, Jl))
                        continue;
                    (void)Jl;
                    add_JTJ_2xN(res, Jp, Jp, 6, 6, s * pix_w());
                }
            }
            else
            {
                for (Index e_idx : pose_odom_edges_[i])
                {
                    const OdomMeas &meas = in_.odom[e_idx];
                    if ((meas.a != i && meas.b != i) || (meas.a != j && meas.b != j))
                        continue;
                    Scalar Ja[6][6], Jb[6][6];
                    odom_jacobian(x, meas, Ja, Jb);
                    if (i == meas.b && j == meas.a)
                        add_JTJ_pose(res, Jb, Ja, s * odom_w());
                    if (i == meas.a && j == meas.b)
                        add_JTJ_pose(res, Ja, Jb, s * odom_w());
                }
            }
        }
        else if (is_ll)
        {
            if (i == j)
            {
                const Index jj = i - n_pose_;
                for (Index e_idx : lmk_obs_edges_[jj])
                {
                    const PixelMeas &obs = in_.observations[e_idx];
                    Pose pose;
                    Landmark lmk;
                    load_pose(x, obs.pose_k, pose);
                    load_lmk(x, obs.lmk_j, lmk);
                    Scalar Jp[2][6], Jl[2][3];
                    if (!project_jacobian(pose, lmk, Jp, Jl))
                        continue;
                    (void)Jp;
                    add_JTJ_2xN(res, Jl, Jl, 3, 3, s * pix_w());
                }
            }
        }
        else if (is_pl || is_lp)
        {
            const Index pose_idx = is_pl ? j : i;
            const Index lmk_idx = is_pl ? (i - n_pose_) : (j - n_pose_);
            for (Index e_idx : lmk_obs_edges_[lmk_idx])
            {
                const PixelMeas &obs = in_.observations[e_idx];
                if (obs.pose_k != pose_idx)
                    continue;
                Pose pose;
                Landmark lmk;
                load_pose(x, obs.pose_k, pose);
                load_lmk(x, obs.lmk_j, lmk);
                Scalar Jp[2][6], Jl[2][3];
                if (!project_jacobian(pose, lmk, Jp, Jl))
                    continue;
                if (is_pl)
                    add_JTJ_2xN(res, Jl, Jp, 3, 6, s * pix_w());
                else
                    add_JTJ_2xN(res, Jp, Jl, 6, 3, s * pix_w());
            }
        }
        return 0;
    }

    // ---- Inspection --------------------------------------------------------
    Index num_observations() const { return static_cast<Index>(in_.observations.size()); }
    Index num_odom() const { return static_cast<Index>(in_.odom.size()); }

    Scalar prior_w() const { return p_.prior_info * cost_scale_; }
    Scalar odom_w() const { return p_.odom_info * cost_scale_; }
    Scalar pix_w() const { return p_.pix_info * cost_scale_; }

private:
    // ---- Block / primal indexing ------------------------------------------
    Index pose_block(Index k) const { return k; }
    Index lmk_block(Index j) const { return n_pose_ + j; }

    // ---- Cost / Jacobian helpers ------------------------------------------

    void load_pose(const Scalar *x, Index k, Pose &p) const
    {
        const Index off = pose_primal_offset(k);
        for (Index a = 0; a < POSE_PRIM_DIM; ++a)
            p[a] = x[off + a];
    }
    void load_lmk(const Scalar *x, Index j, Landmark &l) const
    {
        const Index off = lmk_primal_offset(j);
        for (Index a = 0; a < LMK_DIM; ++a)
            l[a] = x[off + a];
    }

    Scalar prior_sq_residual(const Scalar *x) const
    {
        const Scalar *p0 = x;
        const Scalar *q0 = x + 3;
        Scalar dq[4];
        quat_inv_mul(in_.prior.data() + 3, q0, dq);
        Scalar phi[3];
        quat_log(dq, phi);
        Scalar v = 0.0;
        for (Index a = 0; a < 3; ++a)
            v += (p0[a] - in_.prior[a]) * (p0[a] - in_.prior[a]);
        for (Index a = 0; a < 3; ++a)
            v += phi[a] * phi[a];
        return v;
    }

    void prior_gradient(const Scalar *x, Scalar gp[6]) const
    {
        const Scalar *p0 = x;
        const Scalar *q0 = x + 3;
        for (Index a = 0; a < 3; ++a)
            gp[a] = p0[a] - in_.prior[a];
        Scalar dq[4];
        quat_inv_mul(in_.prior.data() + 3, q0, dq);
        quat_log(dq, gp + 3);
    }

    // Thin wrappers that pull two pose primals out of the global vector
    // and forward to the shared helpers in slam_helpers.hpp.
    void odom_residual(const Scalar *x, const OdomMeas &m, Scalar r[6]) const
    {
        rel_pose_raw(x + pose_primal_offset(m.a), x + pose_primal_offset(m.b), r);
        for (Index a = 0; a < 6; ++a)
            r[a] -= m.z[a];
    }
    void odom_jacobian(const Scalar *x, const OdomMeas &m, Scalar Ja[6][6],
                       Scalar Jb[6][6]) const
    {
        rel_pose_jacobian(x + pose_primal_offset(m.a), x + pose_primal_offset(m.b), Ja, Jb);
    }

    static void add_JTJ_pose(MAT *res, const Scalar A[6][6], const Scalar B[6][6], Scalar s)
    {
        for (Index i = 0; i < 6; ++i)
            for (Index j = 0; j < 6; ++j)
            {
                Scalar v = 0.0;
                for (Index k = 0; k < 6; ++k)
                    v += A[k][i] * B[k][j];
                blasfeo_matel_wrap(res, i, j) += s * v;
            }
    }
    static void add_JT_r_pose(Scalar *res, const Scalar A[6][6], const Scalar r[6], Scalar s)
    {
        for (Index i = 0; i < 6; ++i)
        {
            Scalar v = 0.0;
            for (Index k = 0; k < 6; ++k)
                v += A[k][i] * r[k];
            res[i] += s * v;
        }
    }
    template <typename Acols, typename Bcols>
    static void add_JTJ_2xN(MAT *res, const Acols &A, const Bcols &B, Index m, Index n, Scalar s)
    {
        for (Index i = 0; i < m; ++i)
            for (Index j = 0; j < n; ++j)
            {
                const Scalar v = A[0][i] * B[0][j] + A[1][i] * B[1][j];
                blasfeo_matel_wrap(res, i, j) += s * v;
            }
    }

    void build_edge_indices()
    {
        pose_odom_edges_.assign(n_pose_, {});
        pose_obs_edges_.assign(n_pose_, {});
        lmk_obs_edges_.assign(n_landmark_, {});
        for (Index e = 0; e < static_cast<Index>(in_.odom.size()); ++e)
        {
            pose_odom_edges_[in_.odom[e].a].push_back(e);
            pose_odom_edges_[in_.odom[e].b].push_back(e);
        }
        for (Index e = 0; e < static_cast<Index>(in_.observations.size()); ++e)
        {
            pose_obs_edges_[in_.observations[e].pose_k].push_back(e);
            lmk_obs_edges_[in_.observations[e].lmk_j].push_back(e);
        }
    }

    SlamProblemInput in_;
    SlamProblemParams p_;
    Index n_pose_;
    Index n_landmark_;
    Scalar cost_scale_ = 1.0; // 1 / (total residual count) — see ctor
    std::vector<std::vector<Index>> pose_odom_edges_;
    std::vector<std::vector<Index>> pose_obs_edges_;
    std::vector<std::vector<Index>> lmk_obs_edges_;
};

// ============================================================================
// Driver
// ============================================================================
namespace
{
struct RunResult
{
    IpSolverReturnFlag ret;
    Scalar max_pose_trans_err;
    Scalar max_pose_rot_err;
    Scalar max_lmk_err;
    Scalar elapsed;
};

RunResult solve_and_report(const char *label, bool noisy)
{
    SyntheticSlamParams params;
    params.noisy = noisy;
    SyntheticSlamData data = generate_synthetic_data(params);

    // Hand the problem only what it needs to solve — no ground truth, no
    // knowledge of the synthetic data generator.
    SlamProblemInput input{data.prior, data.odom, data.observations, data.pose_init,
                           data.lmk_init};
    auto problem = std::make_shared<SlamProblem>(std::move(input));
    auto nlp = std::make_shared<NlpGraph>(problem);
    OptionRegistry opts;
    IpAlgBuilder<GraphProblem> builder(nlp);
    auto alg = builder.with_options_registry(&opts).build();
    opts.set_option<Scalar>("tolerance", 1e-6);

    Timer t;
    t.start();
    IpSolverReturnFlag ret = alg->optimize();
    const Scalar elapsed = t.stop();
    auto ipdata = builder.get_ipdata();
    const VecRealView &x = ipdata->current_iterate().primal_x();

    const Index n_pose = data.n_pose();
    Scalar max_pose_trans_err = 0.0;
    Scalar max_pose_rot_err = 0.0;
    for (Index k = 0; k < n_pose; ++k)
    {
        const Index o = problem->pose_primal_offset(k);
        const Scalar dx = x(o + 0) - data.gt_poses[k][0];
        const Scalar dy = x(o + 1) - data.gt_poses[k][1];
        const Scalar dz = x(o + 2) - data.gt_poses[k][2];
        const Scalar t_err = std::sqrt(dx * dx + dy * dy + dz * dz);
        const Scalar q_est[4] = {x(o + 3), x(o + 4), x(o + 5), x(o + 6)};
        Scalar dq[4];
        const Scalar q_gt_inv[4] = {data.gt_poses[k][3], -data.gt_poses[k][4],
                                    -data.gt_poses[k][5], -data.gt_poses[k][6]};
        quat_mul(q_gt_inv, q_est, dq);
        Scalar phi[3];
        quat_log(dq, phi);
        const Scalar r_err = std::sqrt(phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2]);
        max_pose_trans_err = std::max(max_pose_trans_err, t_err);
        max_pose_rot_err = std::max(max_pose_rot_err, r_err);
    }
    const Index n_lmk = data.n_landmark();
    Scalar max_lmk_err = 0.0;
    for (Index j = 0; j < n_lmk; ++j)
    {
        const Index o = problem->lmk_primal_offset(j);
        const Scalar dx = x(o + 0) - data.gt_landmarks[j][0];
        const Scalar dy = x(o + 1) - data.gt_landmarks[j][1];
        const Scalar dz = x(o + 2) - data.gt_landmarks[j][2];
        max_lmk_err = std::max(max_lmk_err, std::sqrt(dx * dx + dy * dy + dz * dz));
    }

    // Loop closures: distinct pose pairs (i, j), separated by at least
    // one full row of the lawn-mower trajectory, that share at least one
    // landmark observation. The trajectory has one row every
    // n_pose / n_rows poses, so two observations of the same landmark
    // across that gap mean the same physical point was seen from two
    // different rows of the survey — exactly the geometric constraint a
    // loop closure provides (and what ties adjacent rows together in
    // aerial-photogrammetry bundle adjustment).
    const Index loop_gap = std::max<Index>(
        1, n_pose / std::max<Index>(1, static_cast<Index>(params.n_rows)));
    std::vector<std::vector<Index>> obs_by_lmk(n_lmk);
    for (const auto &obs : data.observations)
        obs_by_lmk[obs.lmk_j].push_back(obs.pose_k);
    std::set<std::pair<Index, Index>> loop_pairs;
    for (auto &poses_of_lmk : obs_by_lmk)
    {
        std::sort(poses_of_lmk.begin(), poses_of_lmk.end());
        for (std::size_t a = 0; a + 1 < poses_of_lmk.size(); ++a)
            for (std::size_t b = a + 1; b < poses_of_lmk.size(); ++b)
            {
                if (poses_of_lmk[b] - poses_of_lmk[a] >= loop_gap)
                    loop_pairs.emplace(poses_of_lmk[a], poses_of_lmk[b]);
            }
    }

    std::cout << "\n=== fatrop graph SLAM (" << label << ") ===\n";
    std::cout << "  N_pose = " << n_pose << ", N_landmark = " << n_lmk
              << "  (from " << params.n_landmark_candidates << " candidates)\n";
    std::cout << "  measurements: " << problem->num_odom() << " odometry, "
              << problem->num_observations() << " pixel\n";
    std::cout << "  obs / pose: avg = "
              << (problem->num_observations() / static_cast<Scalar>(n_pose)) << "\n";
    std::cout << "  obs / lmk:  avg = "
              << (problem->num_observations() / static_cast<Scalar>(n_lmk)) << "\n";
    std::cout << "  loop closures: " << loop_pairs.size()
              << "  (pose pairs >= " << loop_gap << " apart sharing a landmark)\n";
    std::cout << "  return flag = " << int(ret) << "  elapsed = " << elapsed << " s\n";
    std::cout << "  max pose trans err = " << max_pose_trans_err << " m\n";
    std::cout << "  max pose rot err   = " << max_pose_rot_err << " rad\n";
    std::cout << "  max landmark err   = " << max_lmk_err << " m\n";
    std::cout << "\n--- Timing statistics ---\n";
    std::cout << ipdata->timing_statistics() << "\n";

    return {ret, max_pose_trans_err, max_pose_rot_err, max_lmk_err, elapsed};
}
} // namespace

int main()
{
    const RunResult clean = solve_and_report("noise-free", /*noisy=*/false);
    const RunResult noisy = solve_and_report("noisy measurements", /*noisy=*/true);

    auto converged = [](IpSolverReturnFlag f) {
        return f == IpSolverReturnFlag::Success ||
               f == IpSolverReturnFlag::StopAtAcceptablePoint;
    };

    const Scalar clean_tol = 1e-6;
    const bool clean_ok = converged(clean.ret) && clean.max_pose_trans_err < clean_tol &&
                          clean.max_pose_rot_err < clean_tol && clean.max_lmk_err < clean_tol;

    const Scalar noisy_t_tol = 5e-1;
    const Scalar noisy_r_tol = 1e-1;
    const Scalar noisy_l_tol = 5e-1;
    const bool noisy_ok = converged(noisy.ret) && noisy.max_pose_trans_err < noisy_t_tol &&
                          noisy.max_pose_rot_err < noisy_r_tol && noisy.max_lmk_err < noisy_l_tol;

    std::cout << "\n--- Summary ---\n";
    std::cout << "noise-free recovery: " << (clean_ok ? "PASS" : "FAIL") << "\n";
    std::cout << "noisy recovery:      " << (noisy_ok ? "PASS" : "FAIL") << "\n";
    return (clean_ok && noisy_ok) ? 0 : 1;
}
