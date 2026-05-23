//
// Dense Lie-group example: 3D point-cloud registration via a unit quaternion.
//
// Solved through the @c DenseAbstract / @c NlpDense pipeline: the primal
// variable @c q lives on the unit quaternions (get_nx() == 4) and the search
// direction lives in so(3) (get_nx_tangent() == 3).
//
//     min  L(q) = 0.5 * sum_i || R(q) p_i - q'_i ||^2
//
// over unit quaternions q with the Lie-group retraction
//     retract(q, v) = q (x) Exp_q(v)
// followed by re-normalization. The Gauss-Newton Hessian and gradient w.r.t.
// the body-frame tangent v are derived from
//     r_i(v) = R(q) Exp([v]_x) p_i - q'_i  ~=  r_i(0) - R(q) [p_i]_x v,
// giving J_i = -R(q) [p_i]_x and therefore
//     g = sum_i J_i^T r_i(0) = - sum_i p_i x (R(q)^T q'_i)
//     H = sum_i J_i^T J_i    = sum_i ( |p_i|^2 I - p_i p_i^T ).
//

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include <fatrop/fatrop.hpp>

#include "so3_helpers.hpp"

using namespace fatrop;
using namespace fatrop::examples;

class PointCloudAlignDense : public DenseAbstract
{
public:
    PointCloudAlignDense(std::vector<std::array<Scalar, 3>> source,
                         std::vector<std::array<Scalar, 3>> target,
                         const Scalar (&q_init)[4])
        : source_(std::move(source)), target_(std::move(target))
    {
        for (int i = 0; i < 4; ++i)
            q_init_[i] = q_init[i];
    }

    Index get_nx() const override { return 4; }         // unit quaternion
    Index get_nx_tangent() const override { return 3; } // so(3)
    Index get_ng() const override { return 0; }
    Index get_ng_ineq() const override { return 0; }

    void apply_retraction(const Scalar *q, const Scalar *delta, const Scalar alpha,
                          Scalar *q_next) override
    {
        Scalar v[3] = {alpha * delta[0], alpha * delta[1], alpha * delta[2]};
        Scalar dq[4];
        quat_exp(v, dq);
        Scalar q_new[4];
        quat_mul(q, dq, q_new);
        quat_normalize(q_new);
        for (int i = 0; i < 4; ++i)
            q_next[i] = q_new[i];
    }

    // Gauss-Newton Hessian + gradient packed into the (nx_tan + 1) x nx_tan Hht layout.
    Index eval_Hh(const Scalar *objective_scale, const Scalar *x, const Scalar * /*lam*/,
                  MAT *res) override
    {
        blasfeo_gese_wrap(4, 3, 0.0, res, 0, 0);
        const Scalar q_inv[4] = {x[0], -x[1], -x[2], -x[3]};
        Scalar H[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
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

    Index eval_Ggt(const Scalar *, MAT *res) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        return 0;
    }
    Index eval_Ggt_ineq(const Scalar *, MAT *res) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        return 0;
    }
    Index eval_g(const Scalar *, Scalar *) override { return 0; }
    Index eval_gineq(const Scalar *, Scalar *) override { return 0; }

    Index eval_grad(const Scalar *objective_scale, const Scalar *x, Scalar *res) override
    {
        const Scalar q_inv[4] = {x[0], -x[1], -x[2], -x[3]};
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

    Index eval_f(const Scalar *objective_scale, const Scalar *x, Scalar *res) override
    {
        Scalar sum_sq = 0.;
        for (std::size_t i = 0; i < source_.size(); ++i)
        {
            Scalar Rp[3];
            quat_rotate(x, source_[i].data(), Rp);
            const Scalar dx = Rp[0] - target_[i][0];
            const Scalar dy = Rp[1] - target_[i][1];
            const Scalar dz = Rp[2] - target_[i][2];
            sum_sq += dx * dx + dy * dy + dz * dz;
        }
        *res = 0.5 * objective_scale[0] * sum_sq;
        return 0;
    }

    Index get_bounds(Scalar *, Scalar *) const override { return 0; }
    Index get_initial(Scalar *x) const override
    {
        for (int i = 0; i < 4; ++i)
            x[i] = q_init_[i];
        return 0;
    }

    std::vector<std::array<Scalar, 3>> source_;
    std::vector<std::array<Scalar, 3>> target_;
    Scalar q_init_[4];
};

int main()
{
    std::mt19937 rng(42);
    std::uniform_real_distribution<Scalar> uniform(-1.0, 1.0);
    constexpr Scalar noise_sigma = 0.02;
    std::normal_distribution<Scalar> noise(0.0, noise_sigma);

    constexpr int N = 200;
    std::vector<std::array<Scalar, 3>> cloud_clean(N);
    for (int i = 0; i < N; ++i)
        cloud_clean[i] = {uniform(rng), uniform(rng), uniform(rng)};

    const Scalar axis_raw[3] = {1., 2., 3.};
    const Scalar axis_n = std::sqrt(axis_raw[0] * axis_raw[0] + axis_raw[1] * axis_raw[1] +
                                    axis_raw[2] * axis_raw[2]);
    const Scalar angle = 35.0 * M_PI / 180.0;
    const Scalar half = 0.5 * angle;
    const Scalar sh = std::sin(half) / axis_n;
    const Scalar q_true[4] = {std::cos(half), axis_raw[0] * sh, axis_raw[1] * sh,
                              axis_raw[2] * sh};

    std::vector<std::array<Scalar, 3>> cloud_rotated(N);
    for (int i = 0; i < N; ++i)
    {
        const Scalar pn[3] = {cloud_clean[i][0] + noise(rng), cloud_clean[i][1] + noise(rng),
                              cloud_clean[i][2] + noise(rng)};
        Scalar rp[3];
        quat_rotate(q_true, pn, rp);
        cloud_rotated[i] = {rp[0], rp[1], rp[2]};
    }

    const Scalar q_init[4] = {1., 0., 0., 0.};

    OptionRegistry options;
    auto problem = std::make_shared<PointCloudAlignDense>(cloud_clean, cloud_rotated, q_init);
    IpAlgBuilder<DenseType> builder(std::make_shared<NlpDense>(problem));
    auto ipalg = builder.with_options_registry(&options).build();

    Timer timer;
    timer.start();
    IpSolverReturnFlag ret = ipalg->optimize();
    const auto elapsed = timer.stop();
    auto data = builder.get_ipdata();
    std::cout << "Elapsed time: " << elapsed << "\n";
    std::cout << data->timing_statistics() << "\n";

    const VecRealView &q_final = data->current_iterate().primal_x();
    const Scalar q_est[4] = {q_final(0), q_final(1), q_final(2), q_final(3)};
    std::cout << "Return flag: " << int(ret) << "\n";
    std::cout << "Success:     " << (ret == IpSolverReturnFlag::Success) << "\n";
    std::cout << "q_true  = [" << q_true[0] << ", " << q_true[1] << ", " << q_true[2] << ", "
              << q_true[3] << "]\n";
    std::cout << "q_est   = [" << q_est[0] << ", " << q_est[1] << ", " << q_est[2] << ", "
              << q_est[3] << "]\n";

    const Scalar q_true_inv[4] = {q_true[0], -q_true[1], -q_true[2], -q_true[3]};
    Scalar q_err[4];
    quat_mul(q_true_inv, q_est, q_err);
    Scalar v_err[3];
    quat_log(q_err, v_err);
    const Scalar rot_err =
        std::sqrt(v_err[0] * v_err[0] + v_err[1] * v_err[1] + v_err[2] * v_err[2]);
    std::cout << "rotation error (rad) = " << rot_err << " (~ " << rot_err * 180.0 / M_PI
              << " deg)\n";

    return (ret == IpSolverReturnFlag::Success && rot_err < 0.05) ? 0 : 1;
}
