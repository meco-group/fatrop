//
// Minimal Lie-group OCP example: rotate a unit quaternion onto a target orientation.
//
// The primal state is a unit quaternion in R^4 (Hamilton convention, [w, x, y, z]) and the
// search direction lives in the corresponding Lie algebra so(3) ~= R^3 (the body-frame
// rotation vector). This makes nx = 4 but nx_tangent = 3.
//
// Single-stage OCP (K = 1, nu = 0, no eq/ineq constraints). The cost is the squared
// geodesic distance on SO(3):
//     L(q) = 0.5 * || log( q_target^{-1} (x) q ) ||^2
// with gradient (in the body-frame so(3) at q):
//     grad_v L = log( q_target^{-1} (x) q )
// and a Gauss-Newton Hessian approximation H = I_3.
//
// The retraction uses the standard exponential on the unit quaternions:
//     retract(q, v) = q (x) exp_quat(0.5 * v)
// followed by renormalization for numerical robustness.

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
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

// q_target^{-1} (x) q for unit quaternions: inverse = conjugate.
inline void error_quat(const Scalar *q, const Scalar *q_target, Scalar *q_err)
{
    Scalar q_target_inv[4] = {q_target[0], -q_target[1], -q_target[2], -q_target[3]};
    quat_mul(q_target_inv, q, q_err);
}
} // namespace

class QuaternionLieOcp : public OcpAbstract
{
public:
    QuaternionLieOcp(const Scalar (&q_target)[4], const Scalar (&q_init)[4])
    {
        for (int i = 0; i < 4; ++i)
        {
            q_target_[i] = q_target[i];
            q_init_[i] = q_init[i];
        }
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

    // Gauss-Newton Hessian: I_3 in the tangent space.
    Index eval_RSQrqt(const Scalar *objective_scale, const Scalar * /*inputs_k*/,
                      const Scalar *states_k, const Scalar * /*lam_dyn_k*/,
                      const Scalar * /*lam_eq_k*/, const Scalar * /*lam_eq_ineq_k*/, MAT *res,
                      const Index /*k*/) override
    {
        // Layout: ((nu + nx_tan + 1) x (nu + nx_tan)) = (4 x 3); top 3x3 is the Hessian
        // block, the last row is the gradient.
        blasfeo_gese_wrap(4, 3, 0.0, res, 0, 0);
        blasfeo_diare_wrap(3, objective_scale[0] * 1.0, res, 0, 0);
        // gradient = log(q_target^{-1} (x) q)
        Scalar q_err[4];
        error_quat(states_k, q_target_, q_err);
        Scalar grad[3];
        quat_log(q_err, grad);
        blasfeo_matel_wrap(res, 3, 0) = objective_scale[0] * grad[0];
        blasfeo_matel_wrap(res, 3, 1) = objective_scale[0] * grad[1];
        blasfeo_matel_wrap(res, 3, 2) = objective_scale[0] * grad[2];
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
        Scalar q_err[4];
        error_quat(states_k, q_target_, q_err);
        Scalar grad[3];
        quat_log(q_err, grad);
        res[0] = objective_scale[0] * grad[0];
        res[1] = objective_scale[0] * grad[1];
        res[2] = objective_scale[0] * grad[2];
        return 0;
    }

    Index eval_L(const Scalar *objective_scale, const Scalar * /*inputs_k*/,
                 const Scalar *states_k, Scalar *res, const Index /*k*/) override
    {
        Scalar q_err[4];
        error_quat(states_k, q_target_, q_err);
        Scalar v[3];
        quat_log(q_err, v);
        *res = 0.5 * objective_scale[0] * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
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

    Scalar q_target_[4];
    Scalar q_init_[4];
};

int main()
{
    // Target: 90-degree rotation around z; initial: identity.
    const Scalar q_target[4] = {std::cos(M_PI / 4), 0., 0., std::sin(M_PI / 4)};
    const Scalar q_init[4] = {1., 0., 0., 0.};

    OptionRegistry options;
    auto ocp = std::make_shared<QuaternionLieOcp>(q_target, q_init);
    IpAlgBuilder<OcpType> builder(std::make_shared<NlpOcp>(ocp));
    auto ipalg = builder.with_options_registry(&options).build();

    IpSolverReturnFlag ret = ipalg->optimize();
    auto data = builder.get_ipdata();

    std::cout << "Return flag: " << int(ret) << "\n";
    std::cout << "Success:     " << (ret == IpSolverReturnFlag::Success) << "\n";

    const VecRealView &q_final = data->current_iterate().primal_x();
    std::cout << "q_final  = [" << q_final(0) << ", " << q_final(1) << ", " << q_final(2)
              << ", " << q_final(3) << "]\n";
    std::cout << "q_target = [" << q_target[0] << ", " << q_target[1] << ", " << q_target[2]
              << ", " << q_target[3] << "]\n";

    Scalar q_err[4];
    Scalar q_final_arr[4] = {q_final(0), q_final(1), q_final(2), q_final(3)};
    Scalar q_target_inv[4] = {q_target[0], -q_target[1], -q_target[2], -q_target[3]};
    quat_mul(q_target_inv, q_final_arr, q_err);
    Scalar v[3];
    quat_log(q_err, v);
    const Scalar residual = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    std::cout << "geodesic residual = " << residual << "\n";
    return (ret == IpSolverReturnFlag::Success && residual < 1e-6) ? 0 : 1;
}
