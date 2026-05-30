//
// Quadrotor trajectory optimization on R^3 x R^3 x SO(3) x R^3, with the actual rotor thrusts
// as inputs.  The full model (dynamics, Lie-aware Jacobian / Hessian blocks, cost) lives in the
// shared header quadrotor_lie_model.hpp; see its file header for the derivation.  This file is
// just the driver: a single open-loop trajectory-optimization solve.  For a receding-horizon
// (MPC) loop reusing the same model, see quadrotor_mpc_example.cpp.
//

#include <cmath>
#include <iostream>
#include <memory>

#include "quadrotor_lie_model.hpp"

using namespace fatrop;
using namespace fatrop::examples;

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
