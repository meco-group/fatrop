//
// Mehrotra predictor-corrector QP example.
//
// A 2-D point-mass tracking problem cast as a strictly convex QP and solved
// with the OCP-structured Mehrotra solver from
// fatrop/mehrotra_qp/mehrotra_qp.hxx.
//
// Problem (per stage k = 0..K-1):
//   states  x_k = [px, py, vx, vy]     (4)
//   inputs  u_k = [fx, fy]             (2)
//   dynamics  x_{k+1} = A x_k + B u_k  (linear -- the linearization is exact)
//   stage cost   0.5 * ( q_xy*(px^2 + py^2) + q_v*(vx^2 + vy^2) + r*(fx^2 + fy^2) )
//   terminal cost 0.5 * q_T*( (px - p_ref_x)^2 + (py - p_ref_y)^2 + vx^2 + vy^2 )
//   box constraints on inputs   -u_max <= u_k <= u_max
//   equality at k=0   x_0 = x_init
//
// The objective is quadratic, dynamics and constraints are affine, so this
// is exactly the QP form. We compare the Mehrotra QP solver against the
// fatrop NLP/IP solver as a sanity check.

#include <cmath>
#include <iostream>
#include <memory>

#include <fatrop/fatrop.hpp>

using namespace fatrop;

class TrackingQp : public OcpAbstract
{
public:
    TrackingQp(Index K, Scalar dt, Scalar mass, Scalar u_max, Scalar q_xy, Scalar q_v, Scalar r,
               Scalar q_T, Scalar ref_x, Scalar ref_y)
        : K_(K), dt_(dt), m_(mass), u_max_(u_max), q_xy_(q_xy), q_v_(q_v), r_(r), q_T_(q_T),
          ref_x_(ref_x), ref_y_(ref_y)
    {
    }

    Index get_nx(const Index) const override { return 4; }
    Index get_nu(const Index k) const override { return (k == K_ - 1) ? 0 : 2; }
    // 4 equality constraints at k=0 (initial state pinned), none elsewhere
    Index get_ng(const Index k) const override { return (k == 0) ? 4 : 0; }
    // box constraints on the inputs are inequality constraints
    Index get_ng_ineq(const Index k) const override { return (k == K_ - 1) ? 0 : 2; }
    Index get_horizon_length() const override { return K_; }

    // [B; A; b]^T layout: nu+nx+1 rows, nx columns
    Index eval_BAbt(const Scalar *states_kp1, const Scalar *inputs_k, const Scalar *states_k,
                    MAT *res, const Index k) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        // B (nu=2 rows on top, nx=4 cols):
        //   [  0,    0,      dt/m,   0    ]
        //   [  0,    0,       0,     dt/m ]
        blasfeo_matel_wrap(res, 0, 2) = dt_ / m_;
        blasfeo_matel_wrap(res, 1, 3) = dt_ / m_;
        // A (next 4 rows):
        //   I + dt * [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
        blasfeo_diare_wrap(4, 1.0, res, 2, 0);
        blasfeo_matel_wrap(res, 4, 0) = dt_; // d px / d vx
        blasfeo_matel_wrap(res, 5, 1) = dt_; // d py / d vy
        // b = -x_{k+1} + A x_k + B u_k (last row, used internally)
        Scalar b[4];
        b[0] = -states_kp1[0] + states_k[0] + dt_ * states_k[2];
        b[1] = -states_kp1[1] + states_k[1] + dt_ * states_k[3];
        b[2] = -states_kp1[2] + states_k[2] + dt_ * inputs_k[0] / m_;
        b[3] = -states_kp1[3] + states_k[3] + dt_ * inputs_k[1] / m_;
        for (Index c = 0; c < 4; ++c)
            blasfeo_matel_wrap(res, 6, c) = b[c];
        return 0;
    }

    Index eval_b(const Scalar *states_kp1, const Scalar *inputs_k, const Scalar *states_k,
                 Scalar *res, const Index) override
    {
        res[0] = -states_kp1[0] + states_k[0] + dt_ * states_k[2];
        res[1] = -states_kp1[1] + states_k[1] + dt_ * states_k[3];
        res[2] = -states_kp1[2] + states_k[2] + dt_ * inputs_k[0] / m_;
        res[3] = -states_kp1[3] + states_k[3] + dt_ * inputs_k[1] / m_;
        return 0;
    }

    // Hessian / gradient on top of (u, x). Last row of the matrix is the
    // gradient block; the rest is the Hessian.
    Index eval_RSQrqt(const Scalar *objective_scale, const Scalar *inputs_k,
                      const Scalar *states_k, const Scalar * /*lam_dyn*/,
                      const Scalar * /*lam_eq*/, const Scalar * /*lam_eq_ineq*/, MAT *res,
                      const Index k) override
    {
        const Scalar s = objective_scale[0];
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        if (k == K_ - 1)
        {
            // Hessian: q_T * I_4 on the 4x4 state block
            for (Index i = 0; i < 4; ++i)
                blasfeo_matel_wrap(res, i, i) = s * q_T_;
            // Gradient: q_T*(x - x_target) on the bottom row
            blasfeo_matel_wrap(res, 4, 0) = s * q_T_ * (states_k[0] - ref_x_);
            blasfeo_matel_wrap(res, 4, 1) = s * q_T_ * (states_k[1] - ref_y_);
            blasfeo_matel_wrap(res, 4, 2) = s * q_T_ * states_k[2];
            blasfeo_matel_wrap(res, 4, 3) = s * q_T_ * states_k[3];
        }
        else
        {
            // Hessian = diag(r, r, q_xy, q_xy, q_v, q_v)
            blasfeo_matel_wrap(res, 0, 0) = s * r_;
            blasfeo_matel_wrap(res, 1, 1) = s * r_;
            blasfeo_matel_wrap(res, 2, 2) = s * q_xy_;
            blasfeo_matel_wrap(res, 3, 3) = s * q_xy_;
            blasfeo_matel_wrap(res, 4, 4) = s * q_v_;
            blasfeo_matel_wrap(res, 5, 5) = s * q_v_;
            // Gradient (last row)
            blasfeo_matel_wrap(res, 6, 0) = s * r_ * inputs_k[0];
            blasfeo_matel_wrap(res, 6, 1) = s * r_ * inputs_k[1];
            blasfeo_matel_wrap(res, 6, 2) = s * q_xy_ * states_k[0];
            blasfeo_matel_wrap(res, 6, 3) = s * q_xy_ * states_k[1];
            blasfeo_matel_wrap(res, 6, 4) = s * q_v_ * states_k[2];
            blasfeo_matel_wrap(res, 6, 5) = s * q_v_ * states_k[3];
        }
        return 0;
    }

    Index eval_rq(const Scalar *objective_scale, const Scalar *inputs_k, const Scalar *states_k,
                  Scalar *res, const Index k) override
    {
        const Scalar s = objective_scale[0];
        if (k == K_ - 1)
        {
            res[0] = s * q_T_ * (states_k[0] - ref_x_);
            res[1] = s * q_T_ * (states_k[1] - ref_y_);
            res[2] = s * q_T_ * states_k[2];
            res[3] = s * q_T_ * states_k[3];
        }
        else
        {
            res[0] = s * r_ * inputs_k[0];
            res[1] = s * r_ * inputs_k[1];
            res[2] = s * q_xy_ * states_k[0];
            res[3] = s * q_xy_ * states_k[1];
            res[4] = s * q_v_ * states_k[2];
            res[5] = s * q_v_ * states_k[3];
        }
        return 0;
    }

    Index eval_L(const Scalar *objective_scale, const Scalar *inputs_k, const Scalar *states_k,
                 Scalar *res, const Index k) override
    {
        const Scalar s = objective_scale[0];
        if (k == K_ - 1)
        {
            const Scalar dx = states_k[0] - ref_x_;
            const Scalar dy = states_k[1] - ref_y_;
            *res = 0.5 * s * q_T_ *
                   (dx * dx + dy * dy + states_k[2] * states_k[2] + states_k[3] * states_k[3]);
        }
        else
        {
            *res = 0.5 * s *
                   (r_ * (inputs_k[0] * inputs_k[0] + inputs_k[1] * inputs_k[1]) +
                    q_xy_ * (states_k[0] * states_k[0] + states_k[1] * states_k[1]) +
                    q_v_ * (states_k[2] * states_k[2] + states_k[3] * states_k[3]));
        }
        return 0;
    }

    // Equality jacobian: pins the initial state at k=0
    Index eval_Ggt(const Scalar *, const Scalar *, MAT *res, const Index k) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        if (k == 0)
        {
            // Constraints g_0(x_0) = x_0 - x_init = 0 -> Jacobian w.r.t. x is identity
            blasfeo_diare_wrap(4, 1.0, res, 2, 0); // skip nu=2 rows of input block
        }
        return 0;
    }

    Index eval_g(const Scalar *, const Scalar *states_k, Scalar *res, const Index k) override
    {
        if (k == 0)
        {
            res[0] = states_k[0] - x_init_[0];
            res[1] = states_k[1] - x_init_[1];
            res[2] = states_k[2] - x_init_[2];
            res[3] = states_k[3] - x_init_[3];
        }
        return 0;
    }

    // Inequality jacobian: box constraints on inputs. fatrop expects them as
    //    g_ineq(u, x) in [lower, upper].
    Index eval_Ggt_ineq(const Scalar *, const Scalar *, MAT *res, const Index k) override
    {
        if (k == K_ - 1)
            return 0;
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        blasfeo_matel_wrap(res, 0, 0) = 1.0; // d fx / d u[0]
        blasfeo_matel_wrap(res, 1, 1) = 1.0; // d fy / d u[1]
        return 0;
    }

    Index eval_gineq(const Scalar *inputs_k, const Scalar *, Scalar *res, const Index k) override
    {
        if (k == K_ - 1)
            return 0;
        res[0] = inputs_k[0];
        res[1] = inputs_k[1];
        return 0;
    }

    Index get_bounds(Scalar *lower, Scalar *upper, const Index k) const override
    {
        if (k == K_ - 1)
            return 0;
        lower[0] = -u_max_;
        upper[0] = u_max_;
        lower[1] = -u_max_;
        upper[1] = u_max_;
        return 0;
    }

    Index get_initial_xk(Scalar *xk, const Index) const override
    {
        for (Index i = 0; i < 4; ++i)
            xk[i] = 0.;
        return 0;
    }

    Index get_initial_uk(Scalar *uk, const Index k) const override
    {
        if (k == K_ - 1)
            return 0;
        uk[0] = 0.;
        uk[1] = 0.;
        return 0;
    }

    void set_x_init(Scalar px, Scalar py, Scalar vx, Scalar vy)
    {
        x_init_[0] = px;
        x_init_[1] = py;
        x_init_[2] = vx;
        x_init_[3] = vy;
    }

private:
    Index K_;
    Scalar dt_;
    Scalar m_;
    Scalar u_max_;
    Scalar q_xy_, q_v_, r_, q_T_;
    Scalar ref_x_, ref_y_;
    Scalar x_init_[4] = {0., 0., 0., 0.};
};

int main()
{
    // Build a moderately sized OCP: K=40 stages, dt=0.1s, push the point mass
    // from (0,0) to (2,1) while staying within force bounds.
    const Index K = 40;
    const Scalar dt = 0.1;
    const Scalar mass = 1.0;
    const Scalar u_max = 3.0;
    const Scalar q_xy = 1.0;
    const Scalar q_v = 0.1;
    const Scalar r = 1.0;
    const Scalar q_T = 50.0;
    const Scalar ref_x = 2.0;
    const Scalar ref_y = 1.0;

    auto ocp = std::make_shared<TrackingQp>(K, dt, mass, u_max, q_xy, q_v, r, q_T, ref_x, ref_y);
    ocp->set_x_init(0., 0., 0., 0.);

    // ----- Solve with the Mehrotra predictor-corrector QP solver --------------
    // Usage is intentionally parallel to the NLP solver below: same builder /
    // optimize / get_ipdata / timing_statistics pattern, with
    // MehrotraQpBuilder swapped in for IpAlgBuilder.
    Scalar mehrotra_x_K[4] = {0., 0., 0., 0.};
    Index mehrotra_iters = 0;
    {
        OptionRegistry options;
        MehrotraQpBuilder<OcpType> builder(std::make_shared<NlpOcp>(ocp));
        auto qpalg = builder.with_options_registry(&options).build();
        qpalg->set_tolerance(1e-9);
        qpalg->set_constr_viol_tol(1e-9);
        qpalg->set_max_iter(50);

        Timer timer;
        timer.start();
        IpSolverReturnFlag ret = qpalg->optimize();
        const double elapsed = timer.stop();

        auto data = builder.get_ipdata();
        std::cout << "\n--- Mehrotra QP solver ---\n";
        std::cout << "Elapsed time:    " << elapsed << " s\n";
        std::cout << "Iterations:      " << qpalg->iteration_count() << "\n";
        std::cout << "Final mu:        " << qpalg->final_mu() << "\n";
        std::cout << "Return flag:     " << int(ret) << "  ("
                  << (ret == IpSolverReturnFlag::Success ? "success" : "non-success") << ")\n";

        if (ret != IpSolverReturnFlag::Success)
            return 1;

        const VecRealView &x = data->current_iterate().primal_x();
        const auto &info = data->current_iterate().info();
        const Index off_u0 = info.offsets_primal_u[0];
        const Index off_xN = info.offsets_primal_x[K - 1];
        std::cout << "u_0  = [" << x(off_u0) << ", " << x(off_u0 + 1) << "]\n";
        std::cout << "x_K  = [" << x(off_xN) << ", " << x(off_xN + 1) << ", " << x(off_xN + 2)
                  << ", " << x(off_xN + 3) << "]\n";
        std::cout << "x_ref= [" << ref_x << ", " << ref_y << ", 0, 0]\n";
        for (Index i = 0; i < 4; ++i)
            mehrotra_x_K[i] = x(off_xN + i);
        mehrotra_iters = qpalg->iteration_count();
        std::cout << data->timing_statistics() << std::endl;
    }

    // ----- Same QP solved with the standard fatrop NLP solver for comparison --
    Scalar nlp_x_K[4] = {0., 0., 0., 0.};
    {
        OptionRegistry options;
        IpAlgBuilder<OcpType> builder(std::make_shared<NlpOcp>(ocp));
        auto ipalg = builder.with_options_registry(&options).build();

        Timer timer;
        timer.start();
        IpSolverReturnFlag ret = ipalg->optimize();
        const double elapsed = timer.stop();

        auto data = builder.get_ipdata();
        std::cout << "\n--- fatrop NLP solver (reference) ---\n";
        std::cout << "Elapsed time:    " << elapsed << " s\n";
        std::cout << "Return flag:     " << int(ret) << "  ("
                  << (ret == IpSolverReturnFlag::Success ? "success" : "non-success") << ")\n";
        if (ret != IpSolverReturnFlag::Success)
            return 1;

        const VecRealView &x = data->current_iterate().primal_x();
        const auto &info = data->current_iterate().info();
        const Index off_xN = info.offsets_primal_x[K - 1];
        std::cout << "x_K (NLP) = [" << x(off_xN) << ", " << x(off_xN + 1) << ", "
                  << x(off_xN + 2) << ", " << x(off_xN + 3) << "]\n";
        for (Index i = 0; i < 4; ++i)
            nlp_x_K[i] = x(off_xN + i);
    }

    // Both solvers should land on the same QP optimum (up to convergence
    // tolerance). Compare and report.
    Scalar mismatch = 0.;
    for (Index i = 0; i < 4; ++i)
    {
        const Scalar d = mehrotra_x_K[i] - nlp_x_K[i];
        mismatch = std::max<Scalar>(mismatch, std::abs(d));
    }
    std::cout << "\nMehrotra vs NLP terminal-state max |diff| = " << mismatch
              << "  (Mehrotra needed " << mehrotra_iters << " iterations)\n";

    if (mismatch > 1e-6)
    {
        std::cerr << "Mismatch between Mehrotra and NLP solutions too large.\n";
        return 1;
    }
    std::cout << "Mehrotra QP and fatrop NLP agree -- success.\n";
    return 0;
}
