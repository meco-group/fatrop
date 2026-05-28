//
// Quadrotor trajectory optimization on R^3 x R^3 x SO(3) x R^3, with the actual rotor thrusts
// as inputs.  The four rotor commands f = (f_1, f_2, f_3, f_4) drive the rigid-body motion via
// the standard "+ configuration" mixing
//     T (total thrust) = f_1 + f_2 + f_3 + f_4
//     tau_x            = L (f_2 - f_4)               (roll  -- y-axis rotors)
//     tau_y            = L (f_3 - f_1)               (pitch -- x-axis rotors)
//     tau_z            = c_tau (f_1 - f_2 + f_3 - f_4) (yaw   -- alternating spin)
// so the body angular velocity omega becomes a STATE driven by the rotor-produced torques
// through the Euler equation
//     J omega_dot + omega x J omega = tau(f).
//
// State per stage  x_k = [p (3); v (3); q (4); omega (3)]    nx = 13, nx_tangent = 12
// Input per stage  u_k = f in R^4 (rotor thrusts; no input at terminal)
//
// Continuous-time model:
//     dot p     = v
//     dot v     = R(q) e_3 (T(f) / m) - g e_3
//     dot q     = 0.5 q (x) omega
//     dot omega = J^{-1} ( tau(f) - omega x J omega )
//
// Forward-Euler discretization gives the per-stage defects
//     b_p     = p_k + v_k dt - p_{k+1}
//     b_v     = v_k + (R(q_k) e_3 T(f_k)/m - g e_3) dt - v_{k+1}
//     b_q     = log( q_{k+1}^{-1} q_k (x) Exp(omega_k dt) )
//     b_omega = omega_k + J^{-1} (tau(f_k) - omega_k x J omega_k) dt - omega_{k+1}
//
// The p-, v- and omega-rows are Euclidean (their Jacobian blocks on the next state are already
// -I), so they pass through unchanged. The q-row is the Lie one: its Jacobian on
// delta_theta_{k+1} is -J_l(b_q)^{-1}, not -I. We (implicitly) SCALE that row by M_q = J_l(b_q),
// which puts it in the form fatrop's Riccati recursion expects (coefficient -I on
// delta_theta_{k+1}) with the blocks
//     A_theta = R_{k+1}^T R_k,    B_omega = R_{k+1}^T R_k Exp(omega_k dt) J_r(omega_k dt) dt
// (same as the quaternion example).  The right-hand side is left unchanged by the scaling
// because J_l(b_q) b_q = b_q, so eval_b returns the physical defect and the primal side
// (feasibility, line search, converged trajectory) is exact.
//
// fatrop computes the SCALED dual lambda_tilde of the scaled constraint c_tilde = M c. We use
// lambda_tilde directly in eval_RSQrqt -- no "physical dual recovery" is needed. Justification:
// the row scalings leave the residual unchanged as a FUNCTION (J_l(b_q) b_q = b_q and
// J_r(g_q) g_q = g_q are SO(3) identities, since J_l(v) v = v), so c_tilde and b coincide and
// so do their second derivatives. The doc-recommended "scaled x scaled" contraction
//     sum lambda_tilde_i grad^2 c_tilde_i  =  sum lambda_tilde_i grad^2 b_i
// is therefore exact, and the curvature term we add to the Hessian is just
// sum lambda_tilde_i grad^2 b_i (closed-form for b_v / b_omega, b_q = 0 expansion for b_q).
//
// Closed-form curvature contributions (b_p is linear -> 0):
//   b_v:                   (delta_f, delta_theta_k) cross  +  (delta_theta_k, delta_theta_k)
//   b_omega:               (delta_omega_k, delta_omega_k)
//   b_q  (at b_q -> 0):    (delta_theta_k, delta_omega_k) cross
//                          self-blocks vanish at b_q = 0; O(b_q) corrections dropped.
// The b_q approximation is exact at the optimum, so quadratic convergence is preserved; the
// O(b_q) error is a small relative perturbation away from feasibility.  Initial-condition
// curvature at k=0 vanishes likewise at g_q = 0 and is dropped.
//
// The terminal orientation cost block uses the right-Jacobian-aware GN form (J_r^{-T} J_r^{-1}
// on the theta block, J_r^{-T} dth on the gradient), since
// d/d(delta_theta) Log(q_t^{-1} q Exp(delta_theta))|_0 = J_r(dth)^{-1}.
//
// (We do not add the rotor lower bound f_i >= 0; the cost regularization keeps the thrusts
// near hover, so the relaxed problem matches the physical one for the regimes here.)
//

#include <cmath>
#include <iostream>
#include <memory>
#include <fatrop/fatrop.hpp>

using namespace fatrop;

#include "so3_helpers.hpp"

using namespace fatrop::examples;

namespace
{
// Physical constants for the small-quadrotor instance.
constexpr Scalar G       = 9.81;   // gravity (m / s^2)
constexpr Scalar M_BODY  = 1.0;    // mass (kg)
constexpr Scalar L_ARM   = 0.2;    // rotor arm length (m)
constexpr Scalar C_TAU   = 0.05;   // yaw-torque-per-thrust ratio (m)
constexpr Scalar JX      = 0.005;  // inertia about x_b (kg m^2)
constexpr Scalar JY      = 0.005;  // inertia about y_b
constexpr Scalar JZ      = 0.010;  // inertia about z_b
constexpr Scalar F_HOVER = G * M_BODY / 4.0; // per-rotor hover thrust ( ~2.45 N )

// b_q = log( q_{k+1}^{-1} q_k (x) Exp(omega dt) )  -- physical orientation defect
inline void quat_residual(const Scalar *q_kp1, const Scalar *q_k, const Scalar omega[3],
                          const Scalar dt, Scalar b[3])
{
    Scalar omega_dt[3] = {omega[0] * dt, omega[1] * dt, omega[2] * dt};
    Scalar dq[4]; quat_exp(omega_dt, dq);
    Scalar f[4]; quat_mul(q_k, dq, f);
    Scalar q_inv[4]; quat_inv(q_kp1, q_inv);
    Scalar t[4]; quat_mul(q_inv, f, t);
    quat_log(t, b);
}

// + configuration mixing:  T = sum f_i,  tau = M f
inline void rotor_mixing(const Scalar f[4], Scalar &T, Scalar tau[3])
{
    T      = f[0] + f[1] + f[2] + f[3];
    tau[0] = L_ARM * (f[1] - f[3]);                        // roll
    tau[1] = L_ARM * (f[2] - f[0]);                        // pitch
    tau[2] = C_TAU * (f[0] - f[1] + f[2] - f[3]);          // yaw (alternating spin)
}

} // namespace

class QuadrotorLieOcp : public OcpAbstract
{
public:
    QuadrotorLieOcp(Index K, Scalar dt,
                    const Scalar p_init[3], const Scalar v_init[3], const Scalar q_init[4],
                    const Scalar omega_init[3],
                    const Scalar p_target[3], const Scalar v_target[3], const Scalar q_target[4],
                    Scalar w_f, Scalar w_omega_run,
                    Scalar w_p, Scalar w_v, Scalar w_q, Scalar w_omega_term)
        : K_(K), dt_(dt),
          w_f_(w_f), w_omega_run_(w_omega_run),
          w_p_(w_p), w_v_(w_v), w_q_(w_q), w_omega_term_(w_omega_term)
    {
        for (int i = 0; i < 3; ++i)
        {
            p_init_[i] = p_init[i]; v_init_[i] = v_init[i]; omega_init_[i] = omega_init[i];
            p_target_[i] = p_target[i]; v_target_[i] = v_target[i];
        }
        for (int i = 0; i < 4; ++i) { q_init_[i] = q_init[i]; q_target_[i] = q_target[i]; }
    }

    // ---- dimensions ----
    //   nx = 13 (p3 + v3 + q4 + omega3),  nx_tangent = 12 (delta_p3 + delta_v3 + delta_theta3 + delta_omega3)
    //   nu = 4 (rotor thrusts) on running stages; 0 at terminal.
    //   ng = 12 at k = 0  (fix the initial state: p, v, q, omega).
    Index get_nx(const Index /*k*/) const override { return 13; }
    Index get_nu(const Index k) const override { return k == K_ - 1 ? 0 : 4; }
    Index get_ng(const Index k) const override { return k == 0 ? 12 : 0; }
    Index get_ng_ineq(const Index /*k*/) const override { return 0; }
    Index get_horizon_length() const override { return K_; }
    Index get_nx_tangent(const Index /*k*/) const override { return 12; }
    Index get_nu_tangent(const Index k) const override { return k == K_ - 1 ? 0 : 4; }

    // ---- Retraction: Euclidean on p, v, omega ; right (body-frame) retraction on q ----
    void apply_retraction_xk(const Index /*k*/, const Scalar *xk, const Scalar *delta_xk,
                             const Scalar alpha, Scalar *xk_next) override
    {
        for (int i = 0; i < 3; ++i) xk_next[i]      = xk[i]     + alpha * delta_xk[i];     // p
        for (int i = 0; i < 3; ++i) xk_next[3 + i]  = xk[3 + i] + alpha * delta_xk[3 + i]; // v
        Scalar v[3] = {alpha * delta_xk[6], alpha * delta_xk[7], alpha * delta_xk[8]};     // theta
        Scalar dq[4]; quat_exp(v, dq);
        Scalar qn[4]; quat_mul(xk + 6, dq, qn); quat_normalize(qn);
        for (int i = 0; i < 4; ++i) xk_next[6 + i] = qn[i];
        for (int i = 0; i < 3; ++i) xk_next[10 + i] = xk[10 + i] + alpha * delta_xk[9 + i]; // omega
    }

    // ---- Dynamics linearization ----
    // BAbt layout: (nu_tan + nx_tan + 1) x nx_tan_next = 17 x 12 (input x output, transposed).
    //   rows 0..3   : u_tan = f_1..f_4
    //   rows 4..6   : delta_p_k
    //   rows 7..9   : delta_v_k
    //   rows 10..12 : delta_theta_k
    //   rows 13..15 : delta_omega_k
    //   row  16     : (reserved for b residual)
    //   cols 0..2   : delta_p_{k+1}
    //   cols 3..5   : delta_v_{k+1}
    //   cols 6..8   : delta_theta_{k+1}
    //   cols 9..11  : delta_omega_{k+1}
    //
    // Only the delta_theta_{k+1} block (cols 6..8) is scaled (by J_l(b_q)); the rest are
    // Euclidean and have -I on the next-state side outright.
    Index eval_BAbt(const Scalar *x_kp1, const Scalar *u_k, const Scalar *x_k, MAT *res,
                    const Index /*k*/) override
    {
        blasfeo_gese_wrap(17, 12, 0., res, 0, 0);
        const Scalar *q_k   = x_k   + 6;
        const Scalar *q_kp1 = x_kp1 + 6;
        const Scalar *om_k  = x_k   + 10;
        Scalar R_k[9];   quat_to_rotmat(q_k,   R_k);
        Scalar R_kp1[9]; quat_to_rotmat(q_kp1, R_kp1);
        Scalar T_now; Scalar tau_now[3];
        rotor_mixing(u_k, T_now, tau_now);

        // ---- delta_p_{k+1} block (cols 0..2)  --  depends on delta_p_k (I), delta_v_k (dt I) ----
        blasfeo_diare_wrap(3, 1.0, res, 4, 0);
        blasfeo_diare_wrap(3, dt_, res, 7, 0);

        // ---- delta_v_{k+1} block (cols 3..5) ----
        //   f_i               -> R_k e_3 * dt / m   (each rotor contributes equally via T = sum f_i)
        //   delta_v_k         -> I
        //   delta_theta_k     -> -R_k [e_3]_x T_now dt / m
        const Scalar inv_m_dt = dt_ / M_BODY;
        // f_i rows (0..3) get the same 1x3 row [R_k[r*3+2] * inv_m_dt for r=0..2], so pack a
        // 4x3 column-major block where each column is a constant.
        Scalar f_to_v[12];
        for (int r = 0; r < 3; ++r)
        {
            const Scalar v = R_k[r * 3 + 2] * inv_m_dt;
            f_to_v[r * 4 + 0] = v;
            f_to_v[r * 4 + 1] = v;
            f_to_v[r * 4 + 2] = v;
            f_to_v[r * 4 + 3] = v;
        }
        blasfeo_pack_mat_wrap(4, 3, f_to_v, 4, res, 0, 3);
        blasfeo_diare_wrap(3, 1.0, res, 7, 3);
        Scalar J_v_th[9];
        for (int i = 0; i < 3; ++i)
        {
            J_v_th[i * 3 + 0] =  R_k[i * 3 + 1];   //  R_k * col 1
            J_v_th[i * 3 + 1] = -R_k[i * 3 + 0];   // -R_k * col 0
            J_v_th[i * 3 + 2] = 0.0;
        }
        const Scalar scl_v_th = -T_now * dt_ / M_BODY;
        for (int i = 0; i < 9; ++i) J_v_th[i] *= scl_v_th;
        // J_v_th is stored row-major; reading it column-major (lda = 3) yields its transpose,
        // which is exactly what the BAbt layout (input r, output c) = J_v_th[c, r] needs.
        blasfeo_pack_mat_wrap(3, 3, J_v_th, 3, res, 10, 3);

        // ---- delta_theta_{k+1} block (cols 6..8 ; q-row scaled by J_l(b_q)) ----
        //   A_theta = R_{k+1}^T R_k         on delta_theta_k    (rows 10..12)
        //   B_omega = A_theta Exp(om dt) J_r(om dt) dt  on delta_omega_k (rows 13..15, now a STATE)
        Scalar A_theta[9];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
            {
                Scalar s = 0.;
                for (int k = 0; k < 3; ++k) s += R_kp1[k * 3 + i] * R_k[k * 3 + j];
                A_theta[i * 3 + j] = s;
            }
        Scalar om_dt[3] = {om_k[0] * dt_, om_k[1] * dt_, om_k[2] * dt_};
        Scalar R_dt[9];
        { Scalar om_dt_q[4]; quat_exp(om_dt, om_dt_q); quat_to_rotmat(om_dt_q, R_dt); }
        Scalar Jr_om[9]; so3_right_jacobian(om_dt, Jr_om);
        Scalar tmp[9]; mat3_mul(A_theta, R_dt, tmp);
        Scalar B_om[9]; mat3_mul(tmp, Jr_om, B_om);
        for (int i = 0; i < 9; ++i) B_om[i] *= dt_;
        // Same transpose-by-reinterpretation trick as J_v_th above: pack the row-major buffers
        // as column-major (lda = 3) to land their transpose into BAbt.
        blasfeo_pack_mat_wrap(3, 3, A_theta, 3, res, 10, 6);
        blasfeo_pack_mat_wrap(3, 3, B_om,    3, res, 13, 6);

        // ---- delta_omega_{k+1} block (cols 9..11) ----
        //   f_i           -> dt J^{-1} M[:, i]
        //   delta_omega_k -> I + dt J^{-1} ( [J omega]_x - [omega]_x J )
        const Scalar inv_Jx = 1.0 / JX, inv_Jy = 1.0 / JY, inv_Jz = 1.0 / JZ;
        // Mixing matrix columns (3x4):
        const Scalar Mcol[4][3] = {
            { 0.0,    -L_ARM,  C_TAU },
            { L_ARM,   0.0,   -C_TAU },
            { 0.0,     L_ARM,  C_TAU },
            {-L_ARM,   0.0,   -C_TAU }
        };
        const Scalar inv_J_dt[3] = {inv_Jx * dt_, inv_Jy * dt_, inv_Jz * dt_};
        Scalar f_to_om[12]; // 4x3 column-major: f_to_om[i + 4*j] = (i, j)
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 4; ++i)
                f_to_om[i + 4 * j] = Mcol[i][j] * inv_J_dt[j];
        blasfeo_pack_mat_wrap(4, 3, f_to_om, 4, res, 0, 9);
        // Gyroscopic Jacobian:  M_om = [J omega]_x - [omega]_x J  (with J diagonal).
        const Scalar ox = om_k[0], oy = om_k[1], oz = om_k[2];
        const Scalar M_om_rowmaj[9] = {
            0.0,           (JY - JZ) * oz, (JY - JZ) * oy,
            (JZ - JX) * oz, 0.0,           (JZ - JX) * ox,
            (JX - JY) * oy, (JX - JY) * ox, 0.0
        };
        Scalar A_omega[9];
        for (int j = 0; j < 3; ++j) A_omega[0 * 3 + j] = inv_Jx * M_om_rowmaj[0 * 3 + j] * dt_;
        for (int j = 0; j < 3; ++j) A_omega[1 * 3 + j] = inv_Jy * M_om_rowmaj[1 * 3 + j] * dt_;
        for (int j = 0; j < 3; ++j) A_omega[2 * 3 + j] = inv_Jz * M_om_rowmaj[2 * 3 + j] * dt_;
        A_omega[0] += 1.0; A_omega[4] += 1.0; A_omega[8] += 1.0;
        // A_omega is row-major; reading column-major (lda = 3) writes its transpose, matching
        // the BAbt convention.
        blasfeo_pack_mat_wrap(3, 3, A_omega, 3, res, 13, 9);
        return 0;
    }

    Index eval_b(const Scalar *x_kp1, const Scalar *u_k, const Scalar *x_k, Scalar *res,
                 const Index /*k*/) override
    {
        // Physical defects. The J_l(b_q) row scaling leaves b_q unchanged (J_l(b)b = b), so
        // feasibility / line search see the true defect.
        const Scalar *p_k = x_k,     *v_k = x_k + 3,     *q_k = x_k + 6,     *om_k = x_k + 10;
        const Scalar *p_kp1 = x_kp1, *v_kp1 = x_kp1 + 3, *q_kp1 = x_kp1 + 6, *om_kp1 = x_kp1 + 10;
        Scalar T_now; Scalar tau_now[3];
        rotor_mixing(u_k, T_now, tau_now);

        for (int i = 0; i < 3; ++i) res[i] = p_k[i] + v_k[i] * dt_ - p_kp1[i];
        Scalar R_k[9]; quat_to_rotmat(q_k, R_k);
        for (int i = 0; i < 3; ++i)
            res[3 + i] = v_k[i] + (R_k[i * 3 + 2] * T_now / M_BODY - (i == 2 ? G : 0.0)) * dt_ - v_kp1[i];
        quat_residual(q_kp1, q_k, om_k, dt_, res + 6);
        // b_omega = omega_k + J^{-1}(tau - omega x J omega) dt - omega_{k+1}
        const Scalar Jw[3] = {JX * om_k[0], JY * om_k[1], JZ * om_k[2]};
        const Scalar gyro[3] = {om_k[1] * Jw[2] - om_k[2] * Jw[1],
                                om_k[2] * Jw[0] - om_k[0] * Jw[2],
                                om_k[0] * Jw[1] - om_k[1] * Jw[0]};
        res[9]  = om_k[0] + (tau_now[0] - gyro[0]) / JX * dt_ - om_kp1[0];
        res[10] = om_k[1] + (tau_now[1] - gyro[1]) / JY * dt_ - om_kp1[1];
        res[11] = om_k[2] + (tau_now[2] - gyro[2]) / JZ * dt_ - om_kp1[2];
        return 0;
    }

    // ---- Initial-state path constraint at k=0: (p,v,q,omega) = (init values) ----
    // The q sub-row is scaled by J_r(g_q) (right Jacobian, q_0 un-inverted) so its Jacobian
    // becomes identity.  The residual is unchanged by the scaling (J_r(g) g = g) so eval_g
    // returns the physical g_q and eval_RSQrqt is free to use the scaled dual directly.
    Index eval_Ggt(const Scalar * /*u_k*/, const Scalar * /*x_k*/, MAT *res, const Index k) override
    {
        const Index nu = get_nu_tangent(k);
        const Index nx = get_nx_tangent(k);
        if (k != 0) { blasfeo_gese_wrap(nu + nx + 1, 0, 0., res, 0, 0); return 0; }
        blasfeo_gese_wrap(nu + nx + 1, 12, 0., res, 0, 0);
        // The 12x12 identity sits at row offset nu, splitting across the four state sub-blocks
        // (p, v, theta_scaled, omega). One diare call covers all of them.
        blasfeo_diare_wrap(12, 1.0, res, nu, 0);
        return 0;
    }
    Index eval_g(const Scalar * /*u_k*/, const Scalar *x_k, Scalar *res, const Index k) override
    {
        if (k != 0) return 0;
        for (int i = 0; i < 3; ++i) res[i]     = x_k[i]      - p_init_[i];
        for (int i = 0; i < 3; ++i) res[3 + i] = x_k[3 + i]  - v_init_[i];
        Scalar q_init_inv[4]; quat_inv(q_init_, q_init_inv);
        Scalar t[4]; quat_mul(q_init_inv, x_k + 6, t);
        quat_log(t, res + 6); // J_r(g) g = g, so the physical residual works as the scaled RHS too
        for (int i = 0; i < 3; ++i) res[9 + i] = x_k[10 + i] - omega_init_[i];
        return 0;
    }

    Index eval_Ggt_ineq(const Scalar *, const Scalar *, MAT *res, const Index) override
    { blasfeo_gese_wrap(res->m, res->n, 0., res, 0, 0); return 0; }
    Index eval_gineq(const Scalar *, const Scalar *, Scalar *, const Index) override { return 0; }

    // ---- Cost ----
    // Stage   : 0.5 ( w_f sum (f_i - f_hover)^2 + w_omega_run |omega|^2 )
    // Terminal: 0.5 ( w_p |p - p_target|^2 + w_v |v - v_target|^2
    //               + w_q |Log(q_t^{-1} q)|^2 + w_omega_term |omega|^2 )
    Index eval_L(const Scalar *os, const Scalar *u_k, const Scalar *x_k,
                 Scalar *res, const Index k) override
    {
        if (k == K_ - 1)
        {
            const Scalar *p = x_k, *v = x_k + 3, *q = x_k + 6, *om = x_k + 10;
            const Scalar dp[3] = {p[0] - p_target_[0], p[1] - p_target_[1], p[2] - p_target_[2]};
            const Scalar dv[3] = {v[0] - v_target_[0], v[1] - v_target_[1], v[2] - v_target_[2]};
            Scalar q_t_inv[4]; quat_inv(q_target_, q_t_inv);
            Scalar t[4]; quat_mul(q_t_inv, q, t);
            Scalar dth[3]; quat_log(t, dth);
            *res = 0.5 * os[0] *
                   (w_p_ * (dp[0]*dp[0]   + dp[1]*dp[1]   + dp[2]*dp[2]) +
                    w_v_ * (dv[0]*dv[0]   + dv[1]*dv[1]   + dv[2]*dv[2]) +
                    w_q_ * (dth[0]*dth[0] + dth[1]*dth[1] + dth[2]*dth[2]) +
                    w_omega_term_ * (om[0]*om[0] + om[1]*om[1] + om[2]*om[2]));
        }
        else
        {
            const Scalar *om = x_k + 10;
            Scalar acc = 0.0;
            for (int i = 0; i < 4; ++i) { const Scalar d = u_k[i] - F_HOVER; acc += d * d; }
            *res = 0.5 * os[0] *
                   (w_f_ * acc +
                    w_omega_run_ * (om[0]*om[0] + om[1]*om[1] + om[2]*om[2]));
        }
        return 0;
    }
    Index eval_rq(const Scalar *os, const Scalar *u_k, const Scalar *x_k,
                  Scalar *res, const Index k) override
    {
        if (k == K_ - 1)
        {
            // Gradient w.r.t. x_tan only (12): delta_p, delta_v, delta_theta, delta_omega.
            // For the orientation block: d/d(delta_theta) Log(q_t^{-1} q Exp(delta_theta))|_0
            // = J_r(dth)^{-1}, so d/d(delta_theta) (0.5 w_q |dth|^2) = w_q J_r(dth)^{-T} dth.
            const Scalar *p = x_k, *v = x_k + 3, *q = x_k + 6, *om = x_k + 10;
            for (int i = 0; i < 3; ++i) res[i]     = os[0] * w_p_ * (p[i] - p_target_[i]);
            for (int i = 0; i < 3; ++i) res[3 + i] = os[0] * w_v_ * (v[i] - v_target_[i]);
            Scalar q_t_inv[4]; quat_inv(q_target_, q_t_inv);
            Scalar t[4]; quat_mul(q_t_inv, q, t);
            Scalar dth[3]; quat_log(t, dth);
            Scalar Jr_inv[9]; so3_right_jacobian_inv(dth, Jr_inv);
            for (int i = 0; i < 3; ++i)
            {
                Scalar s = 0.;
                for (int j = 0; j < 3; ++j) s += Jr_inv[j * 3 + i] * dth[j]; // (J_r^{-T} dth)_i
                res[6 + i] = os[0] * w_q_ * s;
            }
            for (int i = 0; i < 3; ++i) res[9 + i] = os[0] * w_omega_term_ * om[i];
        }
        else
        {
            // Gradient: 4 input components, then 12 state-tangent components. Only omega part of
            // the state-tangent is nonzero (running stage cost has no p/v/theta dependence).
            for (int i = 0; i < 4; ++i) res[i] = os[0] * w_f_ * (u_k[i] - F_HOVER);
            for (int i = 0; i < 9; ++i) res[4 + i] = 0.0;
            const Scalar *om = x_k + 10;
            for (int i = 0; i < 3; ++i) res[4 + 9 + i] = os[0] * w_omega_run_ * om[i];
        }
        return 0;
    }

    // Lagrangian Hessian = Gauss-Newton cost block + closed-form within-stage constraint
    // curvature sum lambda_i grad^2 c_i.  We use fatrop's scaled dual lambda_tilde directly
    // (no physical-dual recovery): the row scalings leave the residual unchanged as a
    // function (J_l(b_q) b_q = b_q identically), so c_tilde and b coincide, and likewise
    // their second derivatives -- the doc-recommended "scaled x scaled" contraction
    // Sum lambda_tilde_i grad^2 c_tilde_i is exactly Sum lambda_tilde_i grad^2 b_i.
    // Closed-form pieces (b_p is linear, contributes nothing):
    //   b_v:                  (delta_f, delta_theta) cross + (delta_theta, delta_theta) self
    //   b_omega:              (delta_omega, delta_omega) self (gyroscopic)
    //   b_q  (at b_q = 0):    (delta_theta, delta_omega) cross only
    // b_q self-blocks vanish at b_q = 0; O(b_q) corrections are dropped, which is exact at
    // the optimum (so the tail is still quadratic) and a small relative error elsewhere.
    Index eval_RSQrqt(const Scalar *os, const Scalar *u_k, const Scalar *x_k,
                      const Scalar *lam_dyn_k, const Scalar * /*lam_eq_k*/,
                      const Scalar * /*lam_eq_ineq_k*/, MAT *res, const Index k) override
    {
        const Index nu = get_nu_tangent(k);
        const Index nx = get_nx_tangent(k);
        blasfeo_gese_wrap(nu + nx + 1, nu + nx, 0., res, 0, 0);
        if (k == K_ - 1)
        {
            // Terminal cost block (x-only).  Order: delta_p (0..2), delta_v (3..5),
            // delta_theta (6..8), delta_omega (9..11).  Gradient row at index nx = 12.
            // The orientation block uses the right-Jacobian-aware GN Hessian
            //     H_theta = w_q J_r(dth)^{-T} J_r(dth)^{-1},   grad_theta = w_q J_r(dth)^{-T} dth,
            // which correctly accounts for d/d(delta_theta) Log(q_t^{-1} q Exp(delta_theta))|_0
            // = J_r(dth)^{-1}. At dth = 0 (the optimum) J_r -> I and this reduces to w_q I.
            blasfeo_diare_wrap(3, os[0] * w_p_,          res, 0, 0);
            blasfeo_diare_wrap(3, os[0] * w_v_,          res, 3, 3);
            blasfeo_diare_wrap(3, os[0] * w_omega_term_, res, 9, 9);
            const Scalar *p = x_k, *v = x_k + 3, *q = x_k + 6, *om = x_k + 10;
            Scalar grad_p[3] = {os[0] * w_p_ * (p[0] - p_target_[0]),
                                os[0] * w_p_ * (p[1] - p_target_[1]),
                                os[0] * w_p_ * (p[2] - p_target_[2])};
            Scalar grad_v[3] = {os[0] * w_v_ * (v[0] - v_target_[0]),
                                os[0] * w_v_ * (v[1] - v_target_[1]),
                                os[0] * w_v_ * (v[2] - v_target_[2])};
            blasfeo_pack_mat_wrap(1, 3, grad_p, 1, res, nx, 0);
            blasfeo_pack_mat_wrap(1, 3, grad_v, 1, res, nx, 3);
            Scalar q_t_inv[4]; quat_inv(q_target_, q_t_inv);
            Scalar t[4]; quat_mul(q_t_inv, q, t);
            Scalar dth[3]; quat_log(t, dth);
            Scalar Jr_inv[9]; so3_right_jacobian_inv(dth, Jr_inv);
            // H_theta = w_q (J_r^{-1})^T (J_r^{-1})  -- symmetric 3x3.
            Scalar H_th[9]; // column-major: H_th[i + 3*j] = (i, j)
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                {
                    Scalar s = 0.;
                    for (int kk = 0; kk < 3; ++kk) s += Jr_inv[kk * 3 + i] * Jr_inv[kk * 3 + j];
                    H_th[i + 3 * j] = os[0] * w_q_ * s;
                }
            blasfeo_pack_mat_wrap(3, 3, H_th, 3, res, 6, 6);
            // grad_theta = w_q J_r^{-T} dth
            Scalar grad_th[3];
            for (int i = 0; i < 3; ++i)
            {
                Scalar s = 0.;
                for (int j = 0; j < 3; ++j) s += Jr_inv[j * 3 + i] * dth[j];
                grad_th[i] = os[0] * w_q_ * s;
            }
            blasfeo_pack_mat_wrap(1, 3, grad_th, 1, res, nx, 6);
            Scalar grad_om[3] = {os[0] * w_omega_term_ * om[0],
                                 os[0] * w_omega_term_ * om[1],
                                 os[0] * w_omega_term_ * om[2]};
            blasfeo_pack_mat_wrap(1, 3, grad_om, 1, res, nx, 9);
        }
        else
        {
            // Stage cost block: w_f I on the 4 inputs (0..3); w_omega_run I on the omega
            // sub-block of x_tan (rows/cols nu+9..nu+11).  Gradient row at index nu + nx.
            blasfeo_diare_wrap(4, os[0] * w_f_,         res, 0,        0);
            blasfeo_diare_wrap(3, os[0] * w_omega_run_, res, nu + 9,   nu + 9);
            Scalar grad_f[4] = {os[0] * w_f_ * (u_k[0] - F_HOVER),
                                os[0] * w_f_ * (u_k[1] - F_HOVER),
                                os[0] * w_f_ * (u_k[2] - F_HOVER),
                                os[0] * w_f_ * (u_k[3] - F_HOVER)};
            blasfeo_pack_mat_wrap(1, 4, grad_f, 1, res, nu + nx, 0);
            const Scalar *om = x_k + 10;
            Scalar grad_om[3] = {os[0] * w_omega_run_ * om[0],
                                 os[0] * w_omega_run_ * om[1],
                                 os[0] * w_omega_run_ * om[2]};
            blasfeo_pack_mat_wrap(1, 3, grad_om, 1, res, nu + nx, nu + 9);

            // ---- Closed-form constraint curvature sum lambda_tilde_i grad^2 b_i --------
            // (Uses the scaled dual lambda_tilde directly -- see file header for why this is
            // exact.) b_p is linear -> no contribution. b_v contributes (delta_f, delta_theta)
            // cross and (delta_theta, delta_theta) from the second-order expansion of
            // R(delta_theta) e_3.  b_omega contributes (delta_omega, delta_omega) from the
            // gyroscopic -J^{-1}(omega x J omega) term.  b_q (the Lie-defect Log) contributes a
            // (delta_theta, delta_omega) cross block at b_q = 0:  BCH gives b_q ~ alpha + beta
            // + (1/2)(alpha x beta) with
            //   alpha = R(omega dt)^T delta_theta,   beta = J_r(omega dt) dt delta_omega
            // (forward J_r, since Exp(B + dB) = Exp(B) Exp(J_r(B) dB + O(dB^2))), so
            //   H_{theta, omega}  =  -(dt/2) R(omega dt) [lambda_q]_x J_r(omega dt).
            // The b_q self-blocks vanish at b_q = 0 (eta is linear in dB to O(dB^3), so there
            // is no quadratic-in-delta_omega term to give a self-block).  The O(b_q) corrections
            // we drop are exact at the optimum -- quadratic convergence is preserved -- and a
            // small relative error elsewhere.  Path-constraint curvature at k=0 vanishes
            // likewise at g_q = 0 and is dropped.
            const Scalar *lam_v  = lam_dyn_k + 3;
            const Scalar *lam_q  = lam_dyn_k + 6;
            const Scalar *lam_om = lam_dyn_k + 9;
            const Scalar *q_k_p  = x_k + 6;
            const Scalar *om_k_p = x_k + 10;
            Scalar T_now; Scalar tau_unused[3]; rotor_mixing(u_k, T_now, tau_unused);
            Scalar R_k_p[9]; quat_to_rotmat(q_k_p, R_k_p);
            // w = R(q_k)^T lambda_v
            const Scalar w[3] = {
                R_k_p[0] * lam_v[0] + R_k_p[3] * lam_v[1] + R_k_p[6] * lam_v[2],
                R_k_p[1] * lam_v[0] + R_k_p[4] * lam_v[1] + R_k_p[7] * lam_v[2],
                R_k_p[2] * lam_v[0] + R_k_p[5] * lam_v[1] + R_k_p[8] * lam_v[2]};
            const Scalar inv_m = 1.0 / M_BODY;

            // (delta_f, delta_theta_k) cross from b_v: each f row = (-dt w_1, dt w_0, 0) / m.
            const Scalar fth0 = -dt_ * w[1] * inv_m;
            const Scalar fth1 =  dt_ * w[0] * inv_m;
            for (int j = 0; j < 4; ++j)
            {
                blasfeo_matel_wrap(res, j,        nu + 6) += fth0;
                blasfeo_matel_wrap(res, j,        nu + 7) += fth1;
                blasfeo_matel_wrap(res, nu + 6,   j)      += fth0; // symmetric
                blasfeo_matel_wrap(res, nu + 7,   j)      += fth1;
            }
            // (delta_theta_k, delta_theta_k) from b_v: T*dt/m * [[-w_2, 0, w_0/2],
            //                                                    [0, -w_2, w_1/2],
            //                                                    [w_0/2, w_1/2, 0]]
            const Scalar Tdt_m = T_now * dt_ * inv_m;
            blasfeo_matel_wrap(res, nu + 6, nu + 6) += -w[2] * Tdt_m;
            blasfeo_matel_wrap(res, nu + 7, nu + 7) += -w[2] * Tdt_m;
            const Scalar t_off_02 = 0.5 * w[0] * Tdt_m;
            const Scalar t_off_12 = 0.5 * w[1] * Tdt_m;
            blasfeo_matel_wrap(res, nu + 6, nu + 8) += t_off_02;
            blasfeo_matel_wrap(res, nu + 8, nu + 6) += t_off_02;
            blasfeo_matel_wrap(res, nu + 7, nu + 8) += t_off_12;
            blasfeo_matel_wrap(res, nu + 8, nu + 7) += t_off_12;

            // (delta_omega_k, delta_omega_k) from -dt J^{-1} (omega x J omega).  For diagonal J,
            //   H[a,b] = -dt * (J_b - J_a) * sum_i lambda_om_i * eps_{abi} / J_i
            // (diagonal entries vanish).
            const Scalar mom_01 = -dt_ * (JY - JX) * lam_om[2] / JZ; // = 0 for JX = JY
            const Scalar mom_02 = -dt_ * (JX - JZ) * lam_om[1] / JY;
            const Scalar mom_12 = -dt_ * (JZ - JY) * lam_om[0] / JX;
            blasfeo_matel_wrap(res, nu +  9, nu + 10) += mom_01;
            blasfeo_matel_wrap(res, nu + 10, nu +  9) += mom_01;
            blasfeo_matel_wrap(res, nu +  9, nu + 11) += mom_02;
            blasfeo_matel_wrap(res, nu + 11, nu +  9) += mom_02;
            blasfeo_matel_wrap(res, nu + 10, nu + 11) += mom_12;
            blasfeo_matel_wrap(res, nu + 11, nu + 10) += mom_12;

            // (delta_theta_k, delta_omega_k) cross from b_q (b_q = 0 leading order):
            //   H_{theta, omega} = -(dt/2) R(omega dt) [lambda_q]_x J_r(omega dt)
            // (uses the *forward* right Jacobian, since Exp(B + dB) = Exp(B) Exp(J_r(B) dB + ...)).
            const Scalar om_dt[3] = {om_k_p[0] * dt_, om_k_p[1] * dt_, om_k_p[2] * dt_};
            Scalar R_om[9];
            { Scalar q_om[4]; quat_exp(om_dt, q_om); quat_to_rotmat(q_om, R_om); }
            Scalar lam_q_hat[9]; hat3(lam_q, lam_q_hat);
            Scalar Jr_om[9]; so3_right_jacobian(om_dt, Jr_om);
            Scalar tmp[9]; mat3_mul(R_om, lam_q_hat, tmp);
            Scalar M_qto[9]; mat3_mul(tmp, Jr_om, M_qto);
            for (int i = 0; i < 9; ++i) M_qto[i] *= -0.5 * dt_;
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                {
                    blasfeo_matel_wrap(res, nu + 6 + a, nu + 9 + b) += M_qto[a * 3 + b];
                    blasfeo_matel_wrap(res, nu + 9 + b, nu + 6 + a) += M_qto[a * 3 + b];
                }
        }
        return 0;
    }

    Index get_bounds(Scalar *, Scalar *, const Index) const override { return 0; }
    Index get_initial_xk(Scalar *xk, const Index /*k*/) const override
    {
        for (int i = 0; i < 3; ++i) { xk[i] = p_init_[i]; xk[3 + i] = v_init_[i]; }
        for (int i = 0; i < 4; ++i) xk[6 + i]  = q_init_[i];
        for (int i = 0; i < 3; ++i) xk[10 + i] = omega_init_[i];
        return 0;
    }
    Index get_initial_uk(Scalar *uk, const Index k) const override
    {
        if (k == K_ - 1) return 0;
        for (int i = 0; i < 4; ++i) uk[i] = F_HOVER; // hover thrust per rotor
        return 0;
    }

private:
    const Index K_;
    const Scalar dt_;
    const Scalar w_f_, w_omega_run_;
    const Scalar w_p_, w_v_, w_q_, w_omega_term_;
    Scalar p_init_[3], v_init_[3], q_init_[4], omega_init_[3];
    Scalar p_target_[3], v_target_[3], q_target_[4];
};

int main()
{
    // Translate from origin to (1.5, 1.0, 0.5) m while flipping the body 90 deg about a tilted
    // axis (1,1,0)/sqrt(2). Start and end at rest. Eight 0.25 s stages.
    const Scalar p_init[3]      = {0., 0., 0.};
    const Scalar v_init[3]      = {0., 0., 0.};
    const Scalar q_init[4]      = {1., 0., 0., 0.};
    const Scalar omega_init[3]  = {0., 0., 0.};
    const Scalar p_target[3]    = {1.5, 1.0, 0.5};
    const Scalar v_target[3]    = {0., 0., 0.};
    const Scalar angle          = M_PI / 2.;                                 // 90 deg
    const Scalar ax_n[3]        = {1. / std::sqrt(2.), 1. / std::sqrt(2.), 0.};
    const Scalar sh             = std::sin(angle / 2.);
    const Scalar q_target[4]    = {std::cos(angle / 2.), sh * ax_n[0], sh * ax_n[1], sh * ax_n[2]};
    const Index K  = 8;
    const Scalar dt = 0.25;
    const Scalar w_f = 1.0, w_omega_run = 0.05;
    const Scalar w_p = 50.0, w_v = 5.0, w_q = 50.0, w_omega_term = 5.0;

    auto ocp = std::make_shared<QuadrotorLieOcp>(K, dt, p_init, v_init, q_init, omega_init,
                                                 p_target, v_target, q_target,
                                                 w_f, w_omega_run,
                                                 w_p, w_v, w_q, w_omega_term);
    OptionRegistry options;
    IpAlgBuilder<OcpType> builder(std::make_shared<NlpOcp>(ocp));
    auto ipalg = builder.with_options_registry(&options).build();
    IpSolverReturnFlag ret = ipalg->optimize();
    auto data = builder.get_ipdata();

    const auto &info = data->info();
    const VecRealView &x = data->current_iterate().primal_x();

    const bool converged = (ret == IpSolverReturnFlag::Success ||
                            ret == IpSolverReturnFlag::StopAtAcceptablePoint);
    std::cout << "\n=== Quadrotor Lie OCP (rotor-thrust inputs) ===\n";
    std::cout << "Return flag: " << int(ret) << "  (converged=" << converged << ")\n";
    std::cout << "Iterations:  " << data->iteration_number() << "\n";
    std::cout << "\nTrajectory:\n";
    for (Index k = 0; k < K; ++k)
    {
        const Scalar *xk = x.data() + info.offsets_primal_x[k];
        std::cout << "  k=" << k
                  << "  p=[" << xk[0] << ", " << xk[1] << ", " << xk[2] << "]"
                  << "  v=[" << xk[3] << ", " << xk[4] << ", " << xk[5] << "]"
                  << "  q=[" << xk[6] << ", " << xk[7] << ", " << xk[8] << ", " << xk[9] << "]"
                  << "  omega=[" << xk[10] << ", " << xk[11] << ", " << xk[12] << "]";
        if (k < K - 1)
        {
            const Scalar *uk = x.data() + info.offsets_primal_u[k];
            std::cout << "  f=[" << uk[0] << ", " << uk[1] << ", " << uk[2] << ", " << uk[3] << "]";
        }
        std::cout << "\n";
    }

    const Scalar *x_T = x.data() + info.offsets_primal_x[K - 1];
    const Scalar dp[3] = {x_T[0] - p_target[0], x_T[1] - p_target[1], x_T[2] - p_target[2]};
    Scalar q_t_inv[4]; quat_inv(q_target, q_t_inv);
    Scalar t[4]; quat_mul(q_t_inv, x_T + 6, t);
    Scalar dth[3]; quat_log(t, dth);
    const Scalar dp_norm  = std::sqrt(dp[0]*dp[0]   + dp[1]*dp[1]   + dp[2]*dp[2]);
    const Scalar dth_norm = std::sqrt(dth[0]*dth[0] + dth[1]*dth[1] + dth[2]*dth[2]);
    std::cout << "\nterminal position error    = " << dp_norm  << " m\n";
    std::cout << "terminal orientation error = " << dth_norm << " rad\n";

    // Solver timing breakdown.  compute_function_evaluation() aggregates the time spent inside
    // eval_objective / eval_gradient / eval_constraint_violation / eval_hessian / eval_jacobian
    // (which is where the FD'ed constraint curvature lives for us).
    const auto &timings = data->timing_statistics();
    std::cout << "\n--- timing ---\n" << timings;
    const double tfe = timings.compute_function_evaluation();
    const double tot = timings.full_algorithm.elapsed();
    std::cout << "function evaluation total  = " << tfe << " s   ("
              << (tot > 0. ? 100. * tfe / tot : 0.) << " % of total)\n";
    return (converged && dp_norm < 5e-2 && dth_norm < 5e-2) ? 0 : 1;
}
