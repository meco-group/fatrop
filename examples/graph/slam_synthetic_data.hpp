//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//
// Synthetic large-scale visual-SLAM scenario generation. Pure scenario
// generation — knows nothing about fatrop / NLP / GraphAbstract. Produces:
//
//   * Ground truth: an aerial-drone "lawn-mower" survey trajectory.
//     n_rows parallel rows of poses_per_row poses each are flown at a
//     constant altitude with a fixed downward-looking gimballed camera
//     (the standard photogrammetry mission profile). Adjacent rows are
//     row_spacing m apart in world Y; consecutive poses in a row are
//     step_along_row m apart in world X; every other row is flown in
//     reverse so the trajectory is a single connected chain. Landmarks
//     are scattered on (and slightly above/below) the ground over the
//     survey footprint plus margins. This pattern is realistic for
//     large-scale aerial photogrammetry / mapping and gives a modest
//     number of cross-row loop closures rather than the very-dense
//     covisibility of a tight orbit.
//   * Noisy synthetic measurements: chain odometry between consecutive
//     poses + per-pose pinhole-camera landmark observations (subject to a
//     range / FOV / per-pose-cap filter), with optional Gaussian noise.
//   * Initial guess: dead-reckoned poses (integrating the noisy odometry
//     from the prior) and DLT-triangulated landmarks (using the noisy
//     pixel observations and the dead-reckoned poses — no ground truth
//     dependency).
//
// The SLAM problem layer (graph_slam_example.cpp) consumes the resulting
// SyntheticSlamData by const-reference.
//

#ifndef __fatrop_examples_slam_synthetic_data_hpp__
#define __fatrop_examples_slam_synthetic_data_hpp__

#include "slam_helpers.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <utility>
#include <vector>

namespace fatrop
{
namespace examples
{

// ============================================================================
// Scenario data structures (measurement types live in slam_helpers.hpp)
// ============================================================================
struct SyntheticSlamData
{
    // --- Ground truth ----------------------------------------------------
    std::vector<Pose> gt_poses;          // length params.n_pose
    std::vector<Landmark> gt_landmarks;  // length n_landmark (after filtering)
    Pose prior;                          // anchor prior on pose 0

    // --- Synthetic measurements -----------------------------------------
    std::vector<OdomMeas> odom;          // one per chain edge
    std::vector<PixelMeas> observations; // sparse pose-landmark observations

    // --- Initial guess --------------------------------------------------
    std::vector<Pose> pose_init;     // dead-reckoned poses
    std::vector<Landmark> lmk_init;  // DLT-triangulated landmarks

    Index n_pose() const { return static_cast<Index>(gt_poses.size()); }
    Index n_landmark() const { return static_cast<Index>(gt_landmarks.size()); }
};

struct SyntheticSlamParams
{
    // Aerial drone "lawn-mower" survey at constant altitude. The drone
    // flies n_rows parallel rows of (n_pose / n_rows) poses each; rows
    // are row_spacing m apart in world Y; consecutive poses in a row are
    // step_along_row m apart in world X; every other row is flown in
    // reverse so the trajectory is a single connected chain. The camera
    // is gimballed straight down (z_cam = -Z_world) and shares the same
    // orientation across all poses (typical of photogrammetry rigs).
    Index n_pose = 512;            // total — rounded down to a multiple of n_rows
    Index n_rows = 8;              // → 64 poses per row by default
    Scalar row_spacing = 28.0;     // m between adjacent rows (~27 % sidelap with the default footprint)
    Scalar step_along_row = 6.0;   // m between consecutive poses (~84 % forward overlap)
    Scalar altitude = 12.0;        // m above the ground (constant)

    // Landmark distribution: scattered on (and slightly above/below) the
    // ground over the full survey footprint plus margins so the camera's
    // edge poses still see plenty of landmarks. lmk_z_range models gentle
    // terrain elevation variation.
    Scalar lmk_x_margin = 20.0;
    Scalar lmk_y_margin = 20.0;
    Scalar lmk_z_range = 1.0;
    Scalar view_range = 25.0;       // landmarks farther than this are unseen
    Scalar half_fov_pix = 320.0;    // half-width of the camera image (pixels)
    Index n_landmark_candidates = 3000; // pool size before filtering
    Index max_obs_per_pose = 60;    // cap the per-pose observation count
    Index min_obs_per_lmk = 6;      // drop landmarks observed fewer times

    Scalar odom_t_noise = 3e-3;    // 3 mm per step
    Scalar odom_r_noise = 2e-3;    // 2 mrad per step
    Scalar pix_noise = 0.3;        // 0.3 pixels
    bool noisy = true;

    unsigned rng_seed = 0xC0FFEEu;
};

// ============================================================================
// Internal builders (file-static lambdas — only `generate_synthetic_data` is
// part of the public API).
// ============================================================================
namespace detail
{
// Generate ground-truth trajectory (aerial drone lawn-mower survey at
// constant altitude with a fixed downward-looking gimballed camera) and a
// pool of random 3D landmark candidates uniformly scattered on the ground
// (with mild elevation noise) over the survey footprint plus margins.
inline void gt_world(const SyntheticSlamParams &p, std::mt19937 &rng,
                     std::vector<Pose> &poses, std::vector<Landmark> &lmk_pool)
{
    const Index n_rows = std::max<Index>(1, p.n_rows);
    const Index poses_per_row = std::max<Index>(2, p.n_pose / n_rows);
    const Index n_pose = poses_per_row * n_rows;
    poses.resize(n_pose);
    const Scalar row_length = static_cast<Scalar>(poses_per_row - 1) * p.step_along_row;

    // All poses share a single downward-looking camera orientation
    // (gimballed sensor, typical aerial photogrammetry). z_cam (optical
    // axis) = world -Z; x_cam = world +X; y_cam = world -Y so the image
    // frame is right-handed (x_cam x y_cam = z_cam).
    const Mat3 R_cam{{{1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 0.0, -1.0}}};
    const Quat q_cam = rot_to_quat(R_cam);

    for (Index k = 0; k < n_pose; ++k)
    {
        const Index row = k / poses_per_row;
        const Index in_row = k % poses_per_row;
        const bool reverse = (row % 2 == 1);
        const Scalar t =
            static_cast<Scalar>(in_row) / static_cast<Scalar>(poses_per_row - 1);
        const Scalar x = (reverse ? (1.0 - t) : t) * row_length;
        const Scalar y = static_cast<Scalar>(row) * p.row_spacing;
        poses[k] = {x, y, p.altitude, q_cam[0], q_cam[1], q_cam[2], q_cam[3]};
    }

    // Landmarks on the ground (with a small elevation noise modelling
    // gentle terrain), covering the full survey footprint plus margins.
    const Scalar y_max_row = static_cast<Scalar>(n_rows - 1) * p.row_spacing;
    std::uniform_real_distribution<Scalar> u_x(-p.lmk_x_margin, row_length + p.lmk_x_margin);
    std::uniform_real_distribution<Scalar> u_y(-p.lmk_y_margin, y_max_row + p.lmk_y_margin);
    std::uniform_real_distribution<Scalar> u_z(-p.lmk_z_range, p.lmk_z_range);
    lmk_pool.resize(p.n_landmark_candidates);
    for (Index j = 0; j < p.n_landmark_candidates; ++j)
        lmk_pool[j] = {u_x(rng), u_y(rng), u_z(rng)};
}

// Build chain-odometry measurements (one per consecutive-pose edge),
// optionally noised.
inline std::vector<OdomMeas> make_odometry(const std::vector<Pose> &poses,
                                            const SyntheticSlamParams &p, std::mt19937 &rng)
{
    std::normal_distribution<Scalar> n_t(0.0, p.noisy ? p.odom_t_noise : 0.0);
    std::normal_distribution<Scalar> n_r(0.0, p.noisy ? p.odom_r_noise : 0.0);
    const Index n_pose = static_cast<Index>(poses.size());
    std::vector<OdomMeas> out;
    out.reserve(n_pose - 1);
    for (Index k = 0; k + 1 < n_pose; ++k)
    {
        OdomMeas m;
        m.a = k;
        m.b = k + 1;
        Scalar z[6];
        rel_pose_raw(poses[k].data(), poses[k + 1].data(), z);
        for (Index a = 0; a < 3; ++a)
        {
            m.z[a] = z[a] + n_t(rng);
            m.z[3 + a] = z[3 + a] + n_r(rng);
        }
        out.push_back(m);
    }
    return out;
}

// Build per-pose pinhole observations with FOV / range filters and a
// per-pose cap (closest K in-FOV landmarks). Each observation gets
// pixel-level noise added.
inline std::vector<PixelMeas> make_observations(const std::vector<Pose> &poses,
                                                 const std::vector<Landmark> &lmk_pool,
                                                 const SyntheticSlamParams &p,
                                                 std::mt19937 &rng)
{
    std::normal_distribution<Scalar> n_pix(0.0, p.noisy ? p.pix_noise : 0.0);
    struct Candidate
    {
        Index lmk_j;
        Scalar dist2;
        Scalar u, v;
    };
    std::vector<Candidate> per_pose;
    per_pose.reserve(lmk_pool.size());
    std::vector<PixelMeas> out;
    const Scalar view2 = p.view_range * p.view_range;
    const Index n_pose = static_cast<Index>(poses.size());
    for (Index k = 0; k < n_pose; ++k)
    {
        per_pose.clear();
        const Mat3 Rk = quat_to_rot(poses[k].data() + 3);
        for (Index j = 0; j < static_cast<Index>(lmk_pool.size()); ++j)
        {
            const Vec3 dp{lmk_pool[j][0] - poses[k][0], lmk_pool[j][1] - poses[k][1],
                          lmk_pool[j][2] - poses[k][2]};
            const Scalar d2 = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];
            if (d2 > view2)
                continue;
            const Vec3 lc = rot_t_apply(Rk, dp);
            if (lc[2] <= Z_MIN)
                continue;
            const Scalar u = FOCAL * lc[0] / lc[2];
            const Scalar v = FOCAL * lc[1] / lc[2];
            if (std::abs(u) > p.half_fov_pix || std::abs(v) > p.half_fov_pix)
                continue;
            per_pose.push_back({j, d2, u, v});
        }
        if (static_cast<Index>(per_pose.size()) > p.max_obs_per_pose)
        {
            std::nth_element(per_pose.begin(), per_pose.begin() + p.max_obs_per_pose,
                             per_pose.end(), [](const Candidate &a, const Candidate &b) {
                                 return a.dist2 < b.dist2;
                             });
            per_pose.resize(p.max_obs_per_pose);
        }
        for (const auto &c : per_pose)
        {
            PixelMeas pm;
            pm.pose_k = k;
            pm.lmk_j = c.lmk_j;
            pm.uv[0] = c.u + n_pix(rng);
            pm.uv[1] = c.v + n_pix(rng);
            out.push_back(pm);
        }
    }
    return out;
}

// Drop landmarks with fewer than @c p.min_obs_per_lmk observations, reindex
// the remaining ones and prune the observation list accordingly.
inline void filter_landmarks(SyntheticSlamData &data, const SyntheticSlamParams &p,
                             const std::vector<Landmark> &lmk_pool)
{
    const Index n_cand = static_cast<Index>(lmk_pool.size());
    std::vector<Index> obs_count(n_cand, 0);
    for (const auto &obs : data.observations)
        obs_count[obs.lmk_j]++;
    std::vector<Index> new_idx(n_cand, -1);
    Index next = 0;
    for (Index j = 0; j < n_cand; ++j)
    {
        if (obs_count[j] >= p.min_obs_per_lmk)
        {
            new_idx[j] = next++;
            data.gt_landmarks.push_back(lmk_pool[j]);
        }
    }
    std::vector<PixelMeas> kept;
    kept.reserve(data.observations.size());
    for (const auto &obs : data.observations)
    {
        if (new_idx[obs.lmk_j] < 0)
            continue;
        PixelMeas pm = obs;
        pm.lmk_j = new_idx[obs.lmk_j];
        kept.push_back(pm);
    }
    data.observations = std::move(kept);
}

// Initial guess for poses: integrate the (noisy) odometry from the prior.
// Initial guess for landmarks: linear DLT triangulation using the
// dead-reckoned poses and the (noisy) pixel observations — no
// ground-truth dependency. Each observation gives two linear equations:
//
//      (f * col0(R) - u * col2(R))^T (l - p) = 0
//      (f * col1(R) - v * col2(R))^T (l - p) = 0
//
// Stacking 2*K such equations and solving the normal equations
// (A^T A) l = A^T b with a 3x3 Cholesky gives the multi-view DLT estimate.
inline void build_initial_guess(SyntheticSlamData &data)
{
    const Index n_pose = data.n_pose();
    data.pose_init.resize(n_pose);
    data.pose_init[0] = data.prior;
    for (Index k = 0; k + 1 < n_pose; ++k)
    {
        const Mat3 Rk = quat_to_rot(data.pose_init[k].data() + 3);
        const Scalar *zt = data.odom[k].z;
        const Vec3 dp_world{Rk[0][0] * zt[0] + Rk[0][1] * zt[1] + Rk[0][2] * zt[2],
                            Rk[1][0] * zt[0] + Rk[1][1] * zt[1] + Rk[1][2] * zt[2],
                            Rk[2][0] * zt[0] + Rk[2][1] * zt[1] + Rk[2][2] * zt[2]};
        data.pose_init[k + 1][0] = data.pose_init[k][0] + dp_world[0];
        data.pose_init[k + 1][1] = data.pose_init[k][1] + dp_world[1];
        data.pose_init[k + 1][2] = data.pose_init[k][2] + dp_world[2];
        Scalar dq[4];
        quat_exp(zt + 3, dq);
        Scalar q_new[4];
        quat_mul(data.pose_init[k].data() + 3, dq, q_new);
        quat_normalize(q_new);
        for (Index a = 0; a < 4; ++a)
            data.pose_init[k + 1][3 + a] = q_new[a];
    }

    // Per-landmark observation index.
    const Index n_lmk = data.n_landmark();
    std::vector<std::vector<Index>> by_lmk(n_lmk);
    for (Index e = 0; e < static_cast<Index>(data.observations.size()); ++e)
        by_lmk[data.observations[e].lmk_j].push_back(e);

    data.lmk_init.resize(n_lmk);
    for (Index j = 0; j < n_lmk; ++j)
    {
        Scalar ATA[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        Scalar ATb[3] = {0, 0, 0};
        for (Index e : by_lmk[j])
        {
            const PixelMeas &obs = data.observations[e];
            const Scalar *p_pose = data.pose_init[obs.pose_k].data();
            const Mat3 R = quat_to_rot(p_pose + 3);
            // Columns of R are world-frame directions of camera axes.
            const Vec3 c0{R[0][0], R[1][0], R[2][0]};
            const Vec3 c1{R[0][1], R[1][1], R[2][1]};
            const Vec3 c2{R[0][2], R[1][2], R[2][2]};
            const Vec3 a_u{FOCAL * c0[0] - obs.uv[0] * c2[0],
                           FOCAL * c0[1] - obs.uv[0] * c2[1],
                           FOCAL * c0[2] - obs.uv[0] * c2[2]};
            const Vec3 a_v{FOCAL * c1[0] - obs.uv[1] * c2[0],
                           FOCAL * c1[1] - obs.uv[1] * c2[1],
                           FOCAL * c1[2] - obs.uv[1] * c2[2]};
            const Scalar bu = a_u[0] * p_pose[0] + a_u[1] * p_pose[1] + a_u[2] * p_pose[2];
            const Scalar bv = a_v[0] * p_pose[0] + a_v[1] * p_pose[1] + a_v[2] * p_pose[2];
            for (Index r = 0; r < 3; ++r)
            {
                for (Index c = 0; c < 3; ++c)
                    ATA[r][c] += a_u[r] * a_u[c] + a_v[r] * a_v[c];
                ATb[r] += a_u[r] * bu + a_v[r] * bv;
            }
        }
        // Solve 3x3 SPD system via Cholesky.
        Scalar L[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        bool ok = true;
        for (Index i = 0; i < 3 && ok; ++i)
        {
            Scalar s = ATA[i][i];
            for (Index k = 0; k < i; ++k)
                s -= L[i][k] * L[i][k];
            if (s <= 1e-12)
            {
                ok = false;
                break;
            }
            L[i][i] = std::sqrt(s);
            for (Index rr = i + 1; rr < 3; ++rr)
            {
                Scalar s2 = ATA[rr][i];
                for (Index k = 0; k < i; ++k)
                    s2 -= L[rr][k] * L[i][k];
                L[rr][i] = s2 / L[i][i];
            }
        }
        Landmark l_est{};
        if (ok)
        {
            Scalar y[3];
            for (Index i = 0; i < 3; ++i)
            {
                Scalar s = ATb[i];
                for (Index k = 0; k < i; ++k)
                    s -= L[i][k] * y[k];
                y[i] = s / L[i][i];
            }
            for (Index i = 2; i >= 0; --i)
            {
                Scalar s = y[i];
                for (Index k = i + 1; k < 3; ++k)
                    s -= L[k][i] * l_est[k];
                l_est[i] = s / L[i][i];
            }
        }
        else
        {
            // Degenerate (rank-deficient) case: fall back to the centroid
            // of the observing poses. Rare when min_obs_per_lmk gives
            // landmarks diverse viewpoints.
            for (Index e : by_lmk[j])
                for (Index a = 0; a < 3; ++a)
                    l_est[a] += data.pose_init[data.observations[e].pose_k][a];
            const Scalar w = 1.0 / static_cast<Scalar>(by_lmk[j].size());
            for (Index a = 0; a < 3; ++a)
                l_est[a] *= w;
        }
        data.lmk_init[j] = l_est;
    }
}
} // namespace detail

// ============================================================================
// Public entry point — build a complete synthetic SLAM scenario.
// ============================================================================
inline SyntheticSlamData generate_synthetic_data(const SyntheticSlamParams &p)
{
    std::mt19937 rng_world(p.rng_seed);
    std::mt19937 rng_meas(p.rng_seed + 100);

    SyntheticSlamData data;
    std::vector<Landmark> lmk_pool;
    detail::gt_world(p, rng_world, data.gt_poses, lmk_pool);
    data.prior = data.gt_poses[0];
    data.odom = detail::make_odometry(data.gt_poses, p, rng_meas);
    data.observations = detail::make_observations(data.gt_poses, lmk_pool, p, rng_meas);
    detail::filter_landmarks(data, p, lmk_pool);
    detail::build_initial_guess(data);
    return data;
}

} // namespace examples
} // namespace fatrop

#endif // __fatrop_examples_slam_synthetic_data_hpp__
