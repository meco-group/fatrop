//
// Lie-group OCP example: 3D point-cloud registration via a unit quaternion.
//
// We build a point cloud P = {p_i} in R^3, perturb it with isotropic Gaussian noise to get
// Q = {p_i + eps_i}, rotate the noisy cloud by some unknown rotation R_true to get
// Q' = {R_true (p_i + eps_i)}, and then recover an estimate of R_true by minimizing
//
//     L(q) = 0.5 * sum_i || R(q) p_i - q'_i ||^2
//
// over unit quaternions q. The primal state is q in R^4 (Hamilton convention [w, x, y, z]);
// the search direction lives in the body-frame Lie algebra so(3) ~= R^3. So nx = 4 but
// nx_tangent = 3.
//
// Single-stage OCP (K = 1, nu = 0, no eq/ineq constraints). The Gauss-Newton Hessian and
// gradient w.r.t. the body-frame tangent v are derived from
//     r_i(v) = R(q) exp([v]_x) p_i - q'_i  ~=  r_i(0) - R(q) [p_i]_x v,
// giving J_i = -R(q) [p_i]_x and therefore
//     g = sum_i J_i^T r_i(0)        = - sum_i p_i x (R(q)^T q'_i)
//     H = sum_i J_i^T J_i           = sum_i ( |p_i|^2 I - p_i p_i^T ).
//
// The retraction is the standard exponential on the unit quaternions:
//     retract(q, v) = q (x) exp_quat(0.5 * v),
// followed by renormalization for numerical robustness.

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include <fatrop/fatrop.hpp>

using namespace fatrop;

namespace
{
// Hamilton quaternion product: r = a (x) b, both [w, x, y, z].
inline void quat_mul(const Scalar *a, const Scalar *b, Scalar *r)
{
    const Scalar aw = a[0], ax = a[1], ay = a[2], az = a[3];
    const Scalar bw = b[0], bx = b[1], by = b[2], bz = b[3];
    r[0] = aw * bw - ax * bx - ay * by - az * bz;
    r[1] = aw * bx + ax * bw + ay * bz - az * by;
    r[2] = aw * by - ax * bz + ay * bw + az * bx;
    r[3] = aw * bz + ax * by - ay * bx + az * bw;
}

inline void quat_normalize(Scalar *q)
{
    const Scalar n = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    if (n > 0.)
    {
        q[0] /= n;
        q[1] /= n;
        q[2] /= n;
        q[3] /= n;
    }
}

// exp on the unit quaternions: exp(0.5 * v) where v in R^3 is a rotation vector.
// The factor 1/2 is folded in so callers can pass a tangent step directly.
inline void quat_exp_half(const Scalar *v, Scalar *q)
{
    const Scalar half_angle = 0.5 * std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (half_angle < 1e-12)
    {
        // Second-order Taylor; avoids 0/0.
        q[0] = 1. - 0.5 * half_angle * half_angle;
        q[1] = 0.5 * v[0];
        q[2] = 0.5 * v[1];
        q[3] = 0.5 * v[2];
        quat_normalize(q);
        return;
    }
    const Scalar s = std::sin(half_angle) / (2. * half_angle);
    q[0] = std::cos(half_angle);
    q[1] = s * v[0];
    q[2] = s * v[1];
    q[3] = s * v[2];
}

// log on the unit quaternions: returns the rotation vector v in R^3 such that
// exp_half(v) = q (with sign chosen so the rotation has angle in [-pi, pi]).
inline void quat_log(const Scalar *q, Scalar *v)
{
    Scalar w = q[0];
    Scalar x = q[1], y = q[2], z = q[3];
    // Quaternions q and -q represent the same rotation; pick the half-sphere with w >= 0
    // so the log is well-defined and continuous around the identity.
    if (w < 0.)
    {
        w = -w;
        x = -x;
        y = -y;
        z = -z;
    }
    const Scalar vec_norm = std::sqrt(x * x + y * y + z * z);
    if (vec_norm < 1e-12)
    {
        v[0] = 2. * x;
        v[1] = 2. * y;
        v[2] = 2. * z;
        return;
    }
    const Scalar angle = 2. * std::atan2(vec_norm, w);
    const Scalar s = angle / vec_norm;
    v[0] = s * x;
    v[1] = s * y;
    v[2] = s * z;
}

// Rotate a 3-vector p by a unit quaternion q: p_rot = R(q) p.
// Uses p' = (w^2 - |v|^2) p + 2 (v.p) v + 2 w (v x p) with q = (w, v).
inline void quat_rotate(const Scalar *q, const Scalar *p, Scalar *p_rot)
{
    const Scalar w = q[0], x = q[1], y = q[2], z = q[3];
    const Scalar vp = x * p[0] + y * p[1] + z * p[2];
    const Scalar vv = x * x + y * y + z * z;
    const Scalar c = 1.0 - 2.0 * vv; // = w^2 - |v|^2 for unit q
    const Scalar cx = y * p[2] - z * p[1];
    const Scalar cy = z * p[0] - x * p[2];
    const Scalar cz = x * p[1] - y * p[0];
    p_rot[0] = c * p[0] + 2.0 * vp * x + 2.0 * w * cx;
    p_rot[1] = c * p[1] + 2.0 * vp * y + 2.0 * w * cy;
    p_rot[2] = c * p[2] + 2.0 * vp * z + 2.0 * w * cz;
}

inline void cross3(const Scalar *a, const Scalar *b, Scalar *r)
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}
} // namespace

class PointCloudAlignOcp : public OcpAbstract
{
public:
    PointCloudAlignOcp(std::vector<std::array<Scalar, 3>> source,
                       std::vector<std::array<Scalar, 3>> target,
                       const Scalar (&q_init)[4])
        : source_(std::move(source)), target_(std::move(target))
    {
        for (int i = 0; i < 4; ++i)
            q_init_[i] = q_init[i];
    }

    // Primal (manifold) size = 4 (unit quaternion).
    Index get_nx(const Index /*k*/) const override { return 4; }
    Index get_nu(const Index /*k*/) const override { return 0; }
    Index get_ng(const Index /*k*/) const override { return 0; }
    Index get_ng_ineq(const Index /*k*/) const override { return 0; }
    Index get_horizon_length() const override { return 1; }

    // Tangent (Lie-algebra so(3)) size = 3.
    Index get_nx_tangent(const Index /*k*/) const override { return 3; }
    Index get_nu_tangent(const Index /*k*/) const override { return 0; }

    // Retraction: q_next = q (x) exp(0.5 * alpha * delta) (re-normalized).
    void apply_retraction_xk(const Index /*k*/, const Scalar *xk, const Scalar *delta_xk,
                             const Scalar alpha, Scalar *xk_next) override
    {
        Scalar v[3] = {alpha * delta_xk[0], alpha * delta_xk[1], alpha * delta_xk[2]};
        Scalar dq[4];
        quat_exp_half(v, dq);
        Scalar q_new[4];
        quat_mul(xk, dq, q_new);
        quat_normalize(q_new);
        for (int i = 0; i < 4; ++i)
            xk_next[i] = q_new[i];
    }

    // K = 1, so no dynamics.
    Index eval_BAbt(const Scalar *, const Scalar *, const Scalar *, MAT *, const Index) override
    {
        return 0;
    }

    // Gauss-Newton Hessian H = sum_i (|p_i|^2 I - p_i p_i^T), gradient g = -sum_i p_i x R^T q'_i.
    Index eval_RSQrqt(const Scalar *objective_scale, const Scalar * /*inputs_k*/,
                      const Scalar *states_k, const Scalar * /*lam_dyn_k*/,
                      const Scalar * /*lam_eq_k*/, const Scalar * /*lam_eq_ineq_k*/, MAT *res,
                      const Index /*k*/) override
    {
        // Layout: ((nu + nx_tan + 1) x (nu + nx_tan)) = (4 x 3); top 3x3 is the Hessian
        // block, the last row is the gradient.
        blasfeo_gese_wrap(4, 3, 0.0, res, 0, 0);

        const Scalar q_inv[4] = {states_k[0], -states_k[1], -states_k[2], -states_k[3]};
        Scalar H[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
        Scalar g[3] = {0., 0., 0.};

        for (std::size_t i = 0; i < source_.size(); ++i)
        {
            const Scalar *p = source_[i].data();
            const Scalar *qt = target_[i].data();

            Scalar Rt_qt[3];
            quat_rotate(q_inv, qt, Rt_qt);

            // g -= p_i x (R^T q'_i)
            Scalar cp[3];
            cross3(p, Rt_qt, cp);
            g[0] -= cp[0];
            g[1] -= cp[1];
            g[2] -= cp[2];

            // H += |p_i|^2 I - p_i p_i^T
            const Scalar p2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    H[r][c] += (r == c ? p2 : 0.) - p[r] * p[c];
        }

        const Scalar s = objective_scale[0];
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                blasfeo_matel_wrap(res, r, c) = s * H[r][c];
        for (int c = 0; c < 3; ++c)
            blasfeo_matel_wrap(res, 3, c) = s * g[c];
        return 0;
    }

    // No equality / inequality constraints to linearize.
    Index eval_Ggt(const Scalar *, const Scalar *, MAT *res, const Index) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        return 0;
    }
    Index eval_Ggt_ineq(const Scalar *, const Scalar *, MAT *res, const Index) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        return 0;
    }
    Index eval_b(const Scalar *, const Scalar *, const Scalar *, Scalar *, const Index) override
    {
        return 0;
    }
    Index eval_g(const Scalar *, const Scalar *, Scalar *, const Index) override { return 0; }
    Index eval_gineq(const Scalar *, const Scalar *, Scalar *, const Index) override { return 0; }

    // Objective gradient in tangent space.
    Index eval_rq(const Scalar *objective_scale, const Scalar * /*inputs_k*/,
                  const Scalar *states_k, Scalar *res, const Index /*k*/) override
    {
        const Scalar q_inv[4] = {states_k[0], -states_k[1], -states_k[2], -states_k[3]};
        Scalar g[3] = {0., 0., 0.};
        for (std::size_t i = 0; i < source_.size(); ++i)
        {
            const Scalar *p = source_[i].data();
            const Scalar *qt = target_[i].data();
            Scalar Rt_qt[3];
            quat_rotate(q_inv, qt, Rt_qt);
            Scalar cp[3];
            cross3(p, Rt_qt, cp);
            g[0] -= cp[0];
            g[1] -= cp[1];
            g[2] -= cp[2];
        }
        const Scalar s = objective_scale[0];
        res[0] = s * g[0];
        res[1] = s * g[1];
        res[2] = s * g[2];
        return 0;
    }

    Index eval_L(const Scalar *objective_scale, const Scalar * /*inputs_k*/,
                 const Scalar *states_k, Scalar *res, const Index /*k*/) override
    {
        Scalar sum_sq = 0.;
        for (std::size_t i = 0; i < source_.size(); ++i)
        {
            Scalar Rp[3];
            quat_rotate(states_k, source_[i].data(), Rp);
            const Scalar dx = Rp[0] - target_[i][0];
            const Scalar dy = Rp[1] - target_[i][1];
            const Scalar dz = Rp[2] - target_[i][2];
            sum_sq += dx * dx + dy * dy + dz * dz;
        }
        *res = 0.5 * objective_scale[0] * sum_sq;
        return 0;
    }

    Index get_bounds(Scalar *, Scalar *, const Index) const override { return 0; }

    Index get_initial_xk(Scalar *xk, const Index /*k*/) const override
    {
        for (int i = 0; i < 4; ++i)
            xk[i] = q_init_[i];
        return 0;
    }
    Index get_initial_uk(Scalar *, const Index) const override { return 0; }

    std::vector<std::array<Scalar, 3>> source_;
    std::vector<std::array<Scalar, 3>> target_;
    Scalar q_init_[4];
};

int main()
{
    // ---- Build the point clouds ---------------------------------------------------------
    std::mt19937 rng(42);
    std::uniform_real_distribution<Scalar> uniform(-1.0, 1.0);
    constexpr Scalar noise_sigma = 0.02;
    std::normal_distribution<Scalar> noise(0.0, noise_sigma);

    constexpr int N = 200;
    std::vector<std::array<Scalar, 3>> cloud_clean(N);
    for (int i = 0; i < N; ++i)
        cloud_clean[i] = {uniform(rng), uniform(rng), uniform(rng)};

    // True rotation: 35 degrees around the axis (1, 2, 3)/||.||
    const Scalar axis_raw[3] = {1., 2., 3.};
    const Scalar axis_n = std::sqrt(axis_raw[0] * axis_raw[0] + axis_raw[1] * axis_raw[1] +
                                    axis_raw[2] * axis_raw[2]);
    const Scalar angle = 35.0 * M_PI / 180.0;
    const Scalar half = 0.5 * angle;
    const Scalar sh = std::sin(half) / axis_n;
    const Scalar q_true[4] = {std::cos(half), axis_raw[0] * sh, axis_raw[1] * sh,
                              axis_raw[2] * sh};

    // Second cloud: add isotropic noise to each point, then rotate by q_true.
    std::vector<std::array<Scalar, 3>> cloud_rotated(N);
    for (int i = 0; i < N; ++i)
    {
        const Scalar pn[3] = {cloud_clean[i][0] + noise(rng),
                              cloud_clean[i][1] + noise(rng),
                              cloud_clean[i][2] + noise(rng)};
        Scalar rp[3];
        quat_rotate(q_true, pn, rp);
        cloud_rotated[i] = {rp[0], rp[1], rp[2]};
    }

    // ---- Solve for the rotation that aligns cloud_clean with cloud_rotated --------------
    const Scalar q_init[4] = {1., 0., 0., 0.}; // start from identity

    OptionRegistry options;
    auto ocp = std::make_shared<PointCloudAlignOcp>(cloud_clean, cloud_rotated, q_init);
    IpAlgBuilder<OcpType> builder(std::make_shared<NlpOcp>(ocp));
    auto ipalg = builder.with_options_registry(&options).build();

    Timer timer;
    timer.start();
    IpSolverReturnFlag ret = ipalg->optimize();
    const auto elapsed = timer.stop();
    auto data = builder.get_ipdata();
    std::cout << "Elapsed time: " << elapsed << std::endl;
    std::cout << data->timing_statistics() << std::endl;

    const VecRealView &q_final = data->current_iterate().primal_x();
    const Scalar q_est[4] = {q_final(0), q_final(1), q_final(2), q_final(3)};

    std::cout << "Return flag: " << int(ret) << "\n";
    std::cout << "Success:     " << (ret == IpSolverReturnFlag::Success) << "\n";
    std::cout << "q_true  = [" << q_true[0] << ", " << q_true[1] << ", " << q_true[2] << ", "
              << q_true[3] << "]\n";
    std::cout << "q_est   = [" << q_est[0] << ", " << q_est[1] << ", " << q_est[2] << ", "
              << q_est[3] << "]\n";

    // Geodesic error between estimated and true rotation.
    const Scalar q_true_inv[4] = {q_true[0], -q_true[1], -q_true[2], -q_true[3]};
    Scalar q_err[4];
    quat_mul(q_true_inv, q_est, q_err);
    Scalar v_err[3];
    quat_log(q_err, v_err);
    const Scalar rot_err =
        std::sqrt(v_err[0] * v_err[0] + v_err[1] * v_err[1] + v_err[2] * v_err[2]);
    std::cout << "rotation error (rad) = " << rot_err << " (~ "
              << rot_err * 180.0 / M_PI << " deg)\n";

    // Residual RMSE in point space after applying the recovered rotation.
    Scalar sum_sq = 0.;
    for (int i = 0; i < N; ++i)
    {
        Scalar Rp[3];
        quat_rotate(q_est, cloud_clean[i].data(), Rp);
        const Scalar dx = Rp[0] - cloud_rotated[i][0];
        const Scalar dy = Rp[1] - cloud_rotated[i][1];
        const Scalar dz = Rp[2] - cloud_rotated[i][2];
        sum_sq += dx * dx + dy * dy + dz * dz;
    }
    const Scalar rmse = std::sqrt(sum_sq / N);
    std::cout << "alignment RMSE       = " << rmse << " (noise sigma = " << noise_sigma
              << ")\n";

    // The point-cloud rotation problem is essentially Wahba/Procrustes; with N=200 points
    // and isotropic noise of sigma=0.02 the rotation should recover well within a few
    // hundredths of a radian.
    return (ret == IpSolverReturnFlag::Success && rot_err < 0.05) ? 0 : 1;
}
