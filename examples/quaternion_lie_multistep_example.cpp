//
// Minimal multi-stage Lie-group OCP: rotate a unit-quaternion state from q_init toward a
// target orientation q_target by applying body-frame angular velocities, minimizing
// control effort plus geodesic distance to the target at the final stage.
//
// State (per stage):  q in S^3      unit quaternion [w, x, y, z]
//                     nx = 4 primal, nx_tangent = 3 (so(3))
// Input (k < K-1):    omega in R^3  body-frame angular velocity
//                     nu = 3
//
// Dynamics (forward Euler on SO(3)):
//     q_{k+1} = q_k (x) Exp_q(omega_k * dt)
//
// In so(3) coordinates the dynamics residual is
//     b_k = log( q_{k+1}^{-1} (q_k (x) Exp_q(omega_k * dt)) )^v
// and its body-frame Jacobian on d_theta_{k+1} is -J_l^{-1}(b_k), not -I. To make
// fatrop's Riccati recursion see the standard -I*delta_x_{k+1} + ... form we
// pre-multiply the linearized dynamics row by M = J_l(b_k); the dual hook then
// recovers the physical multiplier J_l(b_k)^T * lambda_tilde.

#include <cmath>
#include <iostream>
#include <memory>
#include <fatrop/fatrop.hpp>

using namespace fatrop;

namespace
{
inline void mat3_identity(Scalar M[9])
{
    for (int i = 0; i < 9; ++i) M[i] = 0.;
    M[0] = M[4] = M[8] = 1.;
}
inline void mat3_mul(const Scalar A[9], const Scalar B[9], Scalar C[9])
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            Scalar s = 0.;
            for (int k = 0; k < 3; ++k) s += A[i * 3 + k] * B[k * 3 + j];
            C[i * 3 + j] = s;
        }
}
inline void hat3(const Scalar v[3], Scalar M[9])
{
    M[0] = 0.;     M[1] = -v[2]; M[2] = v[1];
    M[3] = v[2];   M[4] = 0.;    M[5] = -v[0];
    M[6] = -v[1];  M[7] = v[0];  M[8] = 0.;
}
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
    if (n > 0.) { q[0] /= n; q[1] /= n; q[2] /= n; q[3] /= n; }
}
inline void quat_exp(const Scalar v[3], Scalar q[4])   // q = Exp_q(v) for v in so(3)
{
    const Scalar half = 0.5 * std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (half < 1e-12)
    {
        q[0] = 1. - 0.5 * half * half;
        q[1] = 0.5 * v[0]; q[2] = 0.5 * v[1]; q[3] = 0.5 * v[2];
        quat_normalize(q);
        return;
    }
    const Scalar s = std::sin(half) / (2. * half);
    q[0] = std::cos(half);
    q[1] = s * v[0]; q[2] = s * v[1]; q[3] = s * v[2];
}
inline void quat_log(const Scalar *q, Scalar v[3])
{
    Scalar w = q[0], x = q[1], y = q[2], z = q[3];
    if (w < 0.) { w = -w; x = -x; y = -y; z = -z; }
    const Scalar vec = std::sqrt(x * x + y * y + z * z);
    if (vec < 1e-12) { v[0] = 2. * x; v[1] = 2. * y; v[2] = 2. * z; return; }
    const Scalar angle = 2. * std::atan2(vec, w);
    const Scalar s = angle / vec;
    v[0] = s * x; v[1] = s * y; v[2] = s * z;
}
inline void quat_inv(const Scalar *q, Scalar *qi)
{
    qi[0] = q[0]; qi[1] = -q[1]; qi[2] = -q[2]; qi[3] = -q[3];
}
inline void quat_to_rotmat(const Scalar *q, Scalar R[9])
{
    const Scalar w = q[0], x = q[1], y = q[2], z = q[3];
    R[0] = 1 - 2 * (y * y + z * z); R[1] = 2 * (x * y - z * w); R[2] = 2 * (x * z + y * w);
    R[3] = 2 * (x * y + z * w);     R[4] = 1 - 2 * (x * x + z * z); R[5] = 2 * (y * z - x * w);
    R[6] = 2 * (x * z - y * w);     R[7] = 2 * (y * z + x * w); R[8] = 1 - 2 * (x * x + y * y);
}
// SO(3) left / right Jacobians (closed-form, with the small-angle Taylor near zero).
inline void so3_left_jacobian(const Scalar v[3], Scalar J[9])
{
    const Scalar theta2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    mat3_identity(J);
    Scalar vh[9]; hat3(v, vh);
    Scalar vh2[9]; mat3_mul(vh, vh, vh2);
    if (theta2 < 1e-8)
    {
        for (int i = 0; i < 9; ++i) J[i] += 0.5 * vh[i] + (1. / 6.) * vh2[i];
        return;
    }
    const Scalar theta = std::sqrt(theta2);
    const Scalar a = (1. - std::cos(theta)) / theta2;
    const Scalar b = (theta - std::sin(theta)) / (theta2 * theta);
    for (int i = 0; i < 9; ++i) J[i] += a * vh[i] + b * vh2[i];
}
inline void so3_right_jacobian(const Scalar v[3], Scalar J[9])
{
    const Scalar nv[3] = {-v[0], -v[1], -v[2]};
    so3_left_jacobian(nv, J);
}

// b_k = log( q_{k+1}^{-1} (q_k (x) Exp_q(omega dt)) )^v
inline void compute_residual(const Scalar *q_kp1, const Scalar *q_k, const Scalar omega[3],
                             const Scalar dt, Scalar b[3])
{
    Scalar omega_dt[3] = {omega[0] * dt, omega[1] * dt, omega[2] * dt};
    Scalar dq_omega[4]; quat_exp(omega_dt, dq_omega);
    Scalar f[4]; quat_mul(q_k, dq_omega, f);
    Scalar q_kp1_inv[4]; quat_inv(q_kp1, q_kp1_inv);
    Scalar t[4]; quat_mul(q_kp1_inv, f, t);
    quat_log(t, b);
}
} // namespace

class QuaternionLieMultistepOcp : public OcpAbstract
{
public:
    QuaternionLieMultistepOcp(Index K, Scalar dt, const Scalar q_init[4],
                              const Scalar q_target[4], Scalar omega_weight,
                              Scalar terminal_weight)
        : K_(K), dt_(dt), aO_(omega_weight), bQ_(terminal_weight)
    {
        for (int i = 0; i < 4; ++i)
        {
            q_init_[i] = q_init[i];
            q_target_[i] = q_target[i];
        }
    }

    // ---- dimensions ----
    Index get_nx(const Index /*k*/) const override { return 4; }
    Index get_nu(const Index k) const override { return k == K_ - 1 ? 0 : 3; }
    Index get_ng(const Index k) const override { return k == 0 ? 3 : 0; }
    Index get_ng_ineq(const Index /*k*/) const override { return 0; }
    Index get_horizon_length() const override { return K_; }
    Index get_nx_tangent(const Index /*k*/) const override { return 3; }
    Index get_nu_tangent(const Index k) const override { return k == K_ - 1 ? 0 : 3; }

    // ---- Lie-group retraction on the quaternion state ----
    void apply_retraction_xk(const Index /*k*/, const Scalar *xk, const Scalar *delta_xk,
                             const Scalar alpha, Scalar *xk_next) override
    {
        Scalar v[3] = {alpha * delta_xk[0], alpha * delta_xk[1], alpha * delta_xk[2]};
        Scalar dq[4]; quat_exp(v, dq);
        Scalar q_new[4]; quat_mul(xk, dq, q_new);
        quat_normalize(q_new);
        for (int i = 0; i < 4; ++i) xk_next[i] = q_new[i];
    }

    // ---- Dual transformations: physical_dual = M^T tilde_dual ----
    void apply_dual_eq_dyn_transformation_k(const Index /*k*/, const Scalar *xk,
                                            const Scalar *xkp1, const Scalar *uk,
                                            const Scalar *dual_in, Scalar *dual_out) override
    {
        Scalar b[3]; compute_residual(xkp1, xk, uk, dt_, b);
        Scalar Jl[9]; so3_left_jacobian(b, Jl);
        for (int i = 0; i < 3; ++i)
        {
            Scalar s = 0.;
            for (int j = 0; j < 3; ++j) s += Jl[j * 3 + i] * dual_in[j];
            dual_out[i] = s;
        }
    }
    void apply_dual_eq_path_transformation_k(const Index k, const Scalar * /*uk*/,
                                             const Scalar *xk, const Scalar *dual_in,
                                             Scalar *dual_out) override
    {
        if (k != 0) return;
        Scalar q_init_inv[4]; quat_inv(q_init_, q_init_inv);
        Scalar t[4]; quat_mul(q_init_inv, xk, t);
        Scalar g[3]; quat_log(t, g);
        Scalar Jr[9]; so3_right_jacobian(g, Jr);
        for (int i = 0; i < 3; ++i)
        {
            Scalar s = 0.;
            for (int j = 0; j < 3; ++j) s += Jr[j * 3 + i] * dual_in[j];
            dual_out[i] = s;
        }
    }

    // ---- Dynamics linearization ----
    // BAbt layout: (nu + nx_tan + 1) x nx_tan_next = 7 x 3.
    //   rows 0..2 : d_omega_k        cols 0..2 : b (3 outputs)
    //   rows 3..5 : d_theta_k
    //   row  6    : (reserved for b residual, written by fatrop)
    //
    // After pre-multiplying the row by M = J_l(b), and using the SO(3) identity
    // J_l(b) * J_r^{-1}(b) = Exp(b) on so(3):
    //   A_scaled (d theta_k -> b)   = R_{k+1}^T R_k
    //   B_scaled (d omega_k -> b)   = R_{k+1}^T R_k * Exp(omega_k dt) * J_r(omega_k dt) * dt
    Index eval_BAbt(const Scalar *q_kp1, const Scalar *u_k, const Scalar *q_k, MAT *res,
                    const Index /*k*/) override
    {
        blasfeo_gese_wrap(7, 3, 0., res, 0, 0);
        Scalar R_k[9]; quat_to_rotmat(q_k, R_k);
        Scalar R_kp1[9]; quat_to_rotmat(q_kp1, R_kp1);
        // A_scaled = R_{k+1}^T R_k
        Scalar A[9];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
            {
                Scalar s = 0.;
                for (int k = 0; k < 3; ++k) s += R_kp1[k * 3 + i] * R_k[k * 3 + j];
                A[i * 3 + j] = s;
            }
        // R_dt = Exp(omega dt)  ;  J_r(omega dt)
        Scalar omega_dt[3] = {u_k[0] * dt_, u_k[1] * dt_, u_k[2] * dt_};
        Scalar R_dt[9];
        {
            Scalar omega_dt_q[4]; quat_exp(omega_dt, omega_dt_q);
            quat_to_rotmat(omega_dt_q, R_dt);
        }
        Scalar Jr_omega[9]; so3_right_jacobian(omega_dt, Jr_omega);
        Scalar tmp[9]; mat3_mul(A, R_dt, tmp);
        Scalar B[9]; mat3_mul(tmp, Jr_omega, B);
        for (int i = 0; i < 9; ++i) B[i] *= dt_;
        // write B (rows 0..2) and A (rows 3..5)
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
            {
                blasfeo_matel_wrap(res, 0 + r, c) = B[r * 3 + c];
                blasfeo_matel_wrap(res, 3 + r, c) = A[r * 3 + c];
            }
        return 0;
    }
    Index eval_b(const Scalar *q_kp1, const Scalar *u_k, const Scalar *q_k, Scalar *res,
                 const Index /*k*/) override
    {
        Scalar b[3]; compute_residual(q_kp1, q_k, u_k, dt_, b);
        // Pre-multiply residual by J_l(b) so the linearization matches A_scaled / B_scaled.
        Scalar Jl[9]; so3_left_jacobian(b, Jl);
        for (int i = 0; i < 3; ++i)
        {
            Scalar s = 0.;
            for (int j = 0; j < 3; ++j) s += Jl[i * 3 + j] * b[j];
            res[i] = s;
        }
        return 0;
    }

    // ---- Initial-condition equality: q_0 == q_init (pre-scaled by J_r) ----
    Index eval_Ggt(const Scalar * /*u_k*/, const Scalar * /*x_k*/, MAT *res,
                   const Index k) override
    {
        const Index nu = get_nu_tangent(k);
        const Index nx = get_nx_tangent(k);
        if (k != 0)
        {
            blasfeo_gese_wrap(nu + nx + 1, 0, 0., res, 0, 0);
            return 0;
        }
        blasfeo_gese_wrap(nu + nx + 1, 3, 0., res, 0, 0);
        // After pre-scaling by J_r(g_theta), the Jacobian on d_theta_0 is identity.
        for (int i = 0; i < 3; ++i)
            blasfeo_matel_wrap(res, nu + i, i) = 1.0;
        return 0;
    }
    Index eval_g(const Scalar * /*u_k*/, const Scalar *q_k, Scalar *res, const Index k) override
    {
        if (k != 0) return 0;
        Scalar q_init_inv[4]; quat_inv(q_init_, q_init_inv);
        Scalar t[4]; quat_mul(q_init_inv, q_k, t);
        Scalar g[3]; quat_log(t, g);
        Scalar Jr[9]; so3_right_jacobian(g, Jr);
        for (int i = 0; i < 3; ++i)
        {
            Scalar s = 0.;
            for (int j = 0; j < 3; ++j) s += Jr[i * 3 + j] * g[j];
            res[i] = s;
        }
        return 0;
    }

    Index eval_Ggt_ineq(const Scalar *, const Scalar *, MAT *res, const Index) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0., res, 0, 0);
        return 0;
    }
    Index eval_gineq(const Scalar *, const Scalar *, Scalar *, const Index) override { return 0; }

    // ---- Cost ----
    //   stage k < K-1 : 0.5 * a_omega * ||omega||^2
    //   terminal k = K-1 : 0.5 * b_terminal * ||log(q_target^{-1} q)||^2
    Index eval_L(const Scalar *objective_scale, const Scalar *u_k, const Scalar *q_k,
                 Scalar *res, const Index k) override
    {
        if (k == K_ - 1)
        {
            Scalar q_t_inv[4]; quat_inv(q_target_, q_t_inv);
            Scalar t[4]; quat_mul(q_t_inv, q_k, t);
            Scalar v[3]; quat_log(t, v);
            *res = 0.5 * objective_scale[0] * bQ_ * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        }
        else
        {
            *res = 0.5 * objective_scale[0] * aO_
                   * (u_k[0] * u_k[0] + u_k[1] * u_k[1] + u_k[2] * u_k[2]);
        }
        return 0;
    }
    Index eval_rq(const Scalar *objective_scale, const Scalar *u_k, const Scalar *q_k,
                  Scalar *res, const Index k) override
    {
        if (k == K_ - 1)
        {
            Scalar q_t_inv[4]; quat_inv(q_target_, q_t_inv);
            Scalar t[4]; quat_mul(q_t_inv, q_k, t);
            Scalar v[3]; quat_log(t, v);
            for (int i = 0; i < 3; ++i) res[i] = objective_scale[0] * bQ_ * v[i];
        }
        else
        {
            for (int i = 0; i < 3; ++i) res[i] = objective_scale[0] * aO_ * u_k[i];
            for (int i = 3; i < 6; ++i) res[i] = 0.;
        }
        return 0;
    }
    // Gauss-Newton Hessian: stage cost is a_omega*I on omega; terminal cost is b_terminal*I
    // in the so(3) tangent at the current orientation. No state cost at running stages.
    Index eval_RSQrqt(const Scalar *objective_scale, const Scalar *u_k, const Scalar *q_k,
                      const Scalar * /*lam_dyn_k*/, const Scalar * /*lam_eq_k*/,
                      const Scalar * /*lam_eq_ineq_k*/, MAT *res, const Index k) override
    {
        const Index nu = get_nu_tangent(k);
        const Index nx = get_nx_tangent(k);
        blasfeo_gese_wrap(nu + nx + 1, nu + nx, 0., res, 0, 0);
        if (k == K_ - 1)
        {
            for (int i = 0; i < 3; ++i)
                blasfeo_matel_wrap(res, i, i) = objective_scale[0] * bQ_;
            Scalar q_t_inv[4]; quat_inv(q_target_, q_t_inv);
            Scalar t[4]; quat_mul(q_t_inv, q_k, t);
            Scalar v[3]; quat_log(t, v);
            for (int i = 0; i < 3; ++i)
                blasfeo_matel_wrap(res, nx, i) = objective_scale[0] * bQ_ * v[i];
        }
        else
        {
            for (int i = 0; i < 3; ++i)
                blasfeo_matel_wrap(res, i, i) = objective_scale[0] * aO_;
            for (int i = 0; i < 3; ++i)
                blasfeo_matel_wrap(res, nu + nx, i) = objective_scale[0] * aO_ * u_k[i];
        }
        return 0;
    }

    Index get_bounds(Scalar *, Scalar *, const Index) const override { return 0; }
    Index get_initial_xk(Scalar *xk, const Index /*k*/) const override
    {
        for (int i = 0; i < 4; ++i) xk[i] = q_init_[i];
        return 0;
    }
    Index get_initial_uk(Scalar *uk, const Index k) const override
    {
        if (k == K_ - 1) return 0;
        uk[0] = 0.; uk[1] = 0.; uk[2] = 0.;
        return 0;
    }

private:
    const Index K_;
    const Scalar dt_;
    const Scalar aO_;   // omega regularization weight
    const Scalar bQ_;   // terminal-cost weight on geodesic distance
    Scalar q_init_[4];
    Scalar q_target_[4];
};

int main()
{
    // Rotate from identity to 90 deg around z over 8 stages of 0.25 s each.
    const Scalar q_init[4]   = {1., 0., 0., 0.};
    const Scalar q_target[4] = {std::cos(M_PI / 4.), 0., 0., std::sin(M_PI / 4.)};
    const Index K = 8;
    const Scalar dt = 0.25;
    const Scalar omega_w = 0.1;
    const Scalar terminal_w = 50.0;

    auto ocp = std::make_shared<QuaternionLieMultistepOcp>(K, dt, q_init, q_target,
                                                            omega_w, terminal_w);
    OptionRegistry options;
    IpAlgBuilder<OcpType> builder(std::make_shared<NlpOcp>(ocp));
    auto ipalg = builder.with_options_registry(&options).build();
    IpSolverReturnFlag ret = ipalg->optimize();
    auto data = builder.get_ipdata();

    const auto &info = data->info();
    const VecRealView &x = data->current_iterate().primal_x();

    std::cout << "\n=== Multistep quaternion Lie OCP ===\n";
    std::cout << "Return flag: " << int(ret) << "  (success="
              << (ret == IpSolverReturnFlag::Success) << ")\n";
    std::cout << "\nTrajectory:\n";
    for (Index k = 0; k < K; ++k)
    {
        const Scalar *q = x.data() + info.offsets_primal_x[k];
        std::cout << "  k=" << k << "  q=[" << q[0] << ", " << q[1] << ", " << q[2] << ", "
                  << q[3] << "]";
        if (k < K - 1)
        {
            const Scalar *u = x.data() + info.offsets_primal_u[k];
            std::cout << "  omega=[" << u[0] << ", " << u[1] << ", " << u[2] << "]";
        }
        std::cout << "\n";
    }

    // Final geodesic distance to target.
    const Scalar *q_final = x.data() + info.offsets_primal_x[K - 1];
    Scalar q_t_inv[4]; quat_inv(q_target, q_t_inv);
    Scalar t[4]; quat_mul(q_t_inv, q_final, t);
    Scalar v[3]; quat_log(t, v);
    const Scalar residual = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    std::cout << "\nq_target = [" << q_target[0] << ", " << q_target[1] << ", " << q_target[2]
              << ", " << q_target[3] << "]\n";
    std::cout << "terminal geodesic distance to target = " << residual << "\n";
    return (ret == IpSolverReturnFlag::Success && residual < 1e-2) ? 0 : 1;
}
