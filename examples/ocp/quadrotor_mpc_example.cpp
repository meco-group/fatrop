//
// Closed-loop MPC (receding-horizon) stabilization of the SO(3) quadrotor.
//
// The plant starts from a HEAVY disturbance -- nearly upside-down (a 150 deg flip), spinning
// fast, and thrown off its hover point with a large linear velocity -- and the controller has to
// balance it back to an upright hover at the origin.
//
// Each control step we:
//   1. read the current plant state x (13-dim: p, v, q, omega),
//   2. solve the trajectory-optimization OCP (the exact same QuadrotorLieOcp used by
//      quadrotor_lie_example.cpp) over a short horizon with x as the initial condition and the
//      upright-hover state as the terminal target,
//   3. apply ONLY the first rotor-thrust command f_0 to the plant,
//   4. integrate the true plant forward by one control period (several Euler sub-steps) under
//      that constant f_0,
//   5. shift the horizon by one step and repeat.
//
// The executed closed-loop trajectory is written to quadrotor_mpc.csv for the Python visualizer
// (examples/ocp/quadrotor_mpc_visualize.py).
//
// Each rotor is box-constrained to a physical range  0 <= f_i <= f_max  (an inequality
// constraint in the OCP), so the recovery respects real actuator limits: total thrust is capped
// at a thrust-to-weight ratio of ~2.85 and rotors can never "pull".
//

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

#include "quadrotor_lie_model.hpp"
#include "fatrop/common/printing.hpp"

using namespace fatrop;
using namespace fatrop::examples;

namespace
{
// Orientation error angle (deg) between unit quaternion q and the target q_t.
Scalar tilt_deg(const Scalar q[4], const Scalar q_t[4])
{
    Scalar q_t_inv[4]; quat_inv(q_t, q_t_inv);
    Scalar t[4]; quat_mul(q_t_inv, q, t);
    Scalar dth[3]; quat_log(t, dth);
    return std::sqrt(dth[0]*dth[0] + dth[1]*dth[1] + dth[2]*dth[2]) * 180.0 / M_PI;
}
} // namespace

int main()
{
    // ---- Controller / simulation settings ----
    const Index  K        = 30;     // MPC horizon (stages, incl. terminal)
    const Scalar dt_ctrl  = 0.04;   // control period (s)  -> 1.16 s look-ahead
    const Index  N_steps  = 150;    // number of closed-loop control steps (-> 6.0 s)
    const Index  n_sub    = 8;      // plant integration sub-steps per control period
    const Scalar dt_sub   = dt_ctrl / n_sub;

    // ---- Target: upright hover at the origin, at rest ----
    const Scalar p_tgt[3] = {0., 0., 0.};
    const Scalar v_tgt[3] = {0., 0., 0.};
    const Scalar q_tgt[4] = {1., 0., 0., 0.};
    const Scalar om_tgt[3] = {0., 0., 0.};

    // ---- Heavy initial disturbance ----
    // 150 deg flip about a tilted axis, large body rates, off-target position and velocity.
    Scalar x[13];
    {
        const Scalar p0[3]  = {0.6, -0.4, 0.5};
        const Scalar v0[3]  = {1.2, -1.5, 0.8};
        const Scalar angle  = 150.0 * M_PI / 180.0;
        Scalar ax[3]        = {1.0, 0.4, 0.2};
        const Scalar an     = std::sqrt(ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2]);
        for (int i = 0; i < 3; ++i) ax[i] /= an;
        const Scalar sh     = std::sin(angle / 2.);
        const Scalar q0[4]  = {std::cos(angle / 2.), sh*ax[0], sh*ax[1], sh*ax[2]};
        const Scalar om0[3] = {5.0, -4.0, 3.0};
        for (int i = 0; i < 3; ++i) { x[i] = p0[i]; x[3+i] = v0[i]; x[10+i] = om0[i]; }
        for (int i = 0; i < 4; ++i) x[6+i] = q0[i];
    }

    // ---- Cost weights (running + terminal) ----
    const Scalar w_f = 0.1, w_omega_run = 0.1;
    const Scalar w_p = 100.0, w_v = 20.0, w_q = 200.0, w_omega_term = 25.0;

    // ---- Actuator limits: each rotor f_i in [0, f_max] ----
    // f_max chosen so the four rotors give a max thrust-to-weight ratio of ~3.0.
    const Scalar f_min = 0.0;
    const Scalar f_max = 3.0 * QUAD_G * QUAD_M_BODY / 4.0;

    // Build the OCP once; re-solved every control step with a fresh initial condition.
    auto ocp = std::make_shared<QuadrotorLieOcp>(K, dt_ctrl, x, x + 3, x + 6, x + 10,
                                                 p_tgt, v_tgt, q_tgt,
                                                 w_f, w_omega_run,
                                                 w_p, w_v, w_q, w_omega_term,
                                                 f_min, f_max);

    // Silence the per-iteration solver banner -- the MPC loop solves many OCPs.
    OutputStreamManager::set_stream(std::make_unique<NullStream>());

    std::ofstream csv("quadrotor_mpc.csv");
    csv << std::setprecision(9);
    csv << "step,t,px,py,pz,vx,vy,vz,qw,qx,qy,qz,wx,wy,wz,f0,f1,f2,f3,tilt_deg,omega_norm\n";
    // goal row (step = -1) for the visualizer.
    csv << -1 << ",0,"
        << p_tgt[0] << "," << p_tgt[1] << "," << p_tgt[2] << ",0,0,0,"
        << q_tgt[0] << "," << q_tgt[1] << "," << q_tgt[2] << "," << q_tgt[3]
        << ",0,0,0,,,,,0,0\n";

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "=== Quadrotor SO(3) MPC -- recovery from a heavy disturbance ===\n";
    std::cout << "horizon K=" << K << "  dt=" << dt_ctrl << "s  steps=" << N_steps
              << "  (T=" << N_steps * dt_ctrl << "s)\n";
    std::cout << "initial tilt = " << tilt_deg(x + 6, q_tgt) << " deg,  |omega| = "
              << std::sqrt(x[10]*x[10] + x[11]*x[11] + x[12]*x[12]) << " rad/s\n\n";
    std::cout << " step    t     |p|     |v|    tilt(deg)  |omega|   iters\n";

    Index total_iters = 0;
    Index worst_iters = 0;
    // Solver wall-time accounting across the whole closed-loop run.
    double t_solve_total = 0., t_solve_worst = 0.;
    double t_feval_total = 0., t_sd_total = 0.;
    for (Index step = 0; step < N_steps; ++step)
    {
        const Scalar t_now = step * dt_ctrl;

        // ---- Solve the OCP for the current state ----
        ocp->set_init(x, x + 3, x + 6, x + 10);
        OptionRegistry options;
        IpAlgBuilder<OcpType> builder(std::make_shared<NlpOcp>(ocp));
        auto ipalg = builder.with_options_registry(&options).build();
        // Safety cap: with warm starting solves take a handful of iterations; the cap just keeps
        // a pathological instance from running away.
        options.set_option("max_iter", Index(80));
        IpSolverReturnFlag ret = ipalg->optimize();
        auto data = builder.get_ipdata();
        const Index iters = data->iteration_number();
        total_iters += iters;
        worst_iters = std::max(worst_iters, iters);

        // ---- Per-solve timing (this OCP's own internal timers) ----
        const auto &tm = data->timing_statistics();
        const double t_solve = tm.full_algorithm.elapsed();
        t_solve_total += t_solve;
        t_solve_worst = std::max(t_solve_worst, t_solve);
        t_feval_total += tm.compute_function_evaluation();
        t_sd_total    += tm.compute_search_dir.elapsed();

        const auto &info = data->info();
        const VecRealView &xsol = data->current_iterate().primal_x();
        const Scalar *u0 = xsol.data() + info.offsets_primal_u[0];
        Scalar f[4] = {u0[0], u0[1], u0[2], u0[3]};
        if (ret != IpSolverReturnFlag::Success && ret != IpSolverReturnFlag::StopAtAcceptablePoint)
        {
            // Fall back to hover thrust if a solve hiccups; keeps the loop running.
            for (int i = 0; i < 4; ++i) f[i] = QUAD_F_HOVER;
        }

        // ---- Log the current (pre-step) state and the applied command ----
        const Scalar pn = std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        const Scalar vn = std::sqrt(x[3]*x[3] + x[4]*x[4] + x[5]*x[5]);
        const Scalar on = std::sqrt(x[10]*x[10] + x[11]*x[11] + x[12]*x[12]);
        const Scalar td = tilt_deg(x + 6, q_tgt);
        csv << step << "," << t_now << ",";
        for (int i = 0; i < 13; ++i) csv << x[i] << ",";
        for (int i = 0; i < 4; ++i)  csv << f[i] << ",";
        csv << td << "," << on << "\n";

        if (step % 5 == 0 || step == N_steps - 1)
            std::cout << std::setw(5) << step << " " << std::setw(6) << t_now << " "
                      << std::setw(7) << pn << " " << std::setw(7) << vn << " "
                      << std::setw(9) << td << " " << std::setw(8) << on << " "
                      << std::setw(7) << iters << "\n";

        // ---- Warm-start the next solve with this solution shifted one stage ----
        // (states x[1..K-1] -> x[0..K-2], duplicate the terminal; same for the K-1 controls).
        for (Index k = 0; k + 1 < K; ++k)
            ocp->set_warm_xk(k, xsol.data() + info.offsets_primal_x[k + 1]);
        ocp->set_warm_xk(K - 1, xsol.data() + info.offsets_primal_x[K - 1]);
        for (Index k = 0; k + 2 < K; ++k)
            ocp->set_warm_uk(k, xsol.data() + info.offsets_primal_u[k + 1]);
        ocp->set_warm_uk(K - 2, xsol.data() + info.offsets_primal_u[K - 2]);
        ocp->enable_warm_start(true);

        // ---- Advance the true plant by one control period under constant f ----
        for (Index s = 0; s < n_sub; ++s)
        {
            Scalar x_next[13];
            quad_simulate_step(x, f, dt_sub, x_next);
            for (int i = 0; i < 13; ++i) x[i] = x_next[i];
        }
    }

    // Final state row (step = N_steps) so the visualizer sees the settled pose.
    {
        const Scalar td = tilt_deg(x + 6, q_tgt);
        const Scalar on = std::sqrt(x[10]*x[10] + x[11]*x[11] + x[12]*x[12]);
        csv << N_steps << "," << N_steps * dt_ctrl << ",";
        for (int i = 0; i < 13; ++i) csv << x[i] << ",";
        csv << QUAD_F_HOVER << "," << QUAD_F_HOVER << "," << QUAD_F_HOVER << "," << QUAD_F_HOVER
            << "," << td << "," << on << "\n";
    }
    csv.close();

    const Scalar final_p = std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    const Scalar final_tilt = tilt_deg(x + 6, q_tgt);
    const Scalar final_om = std::sqrt(x[10]*x[10] + x[11]*x[11] + x[12]*x[12]);
    std::cout << "\nfinal: |p| = " << final_p << " m,  tilt = " << final_tilt
              << " deg,  |omega| = " << final_om << " rad/s\n";
    std::cout << "solver: avg " << (double)total_iters / N_steps << " iters/step, worst "
              << worst_iters << "\n";

    // ---- Solver timing report (sum of each step's internal full_algorithm timer) ----
    const double ms = 1e3;
    std::cout << "\n--- solver timings (over " << N_steps << " MPC solves) ---\n"
              << std::setprecision(3)
              << "  total solve wall time : " << t_solve_total * ms << " ms\n"
              << "  avg  per control step : " << t_solve_total / N_steps * ms << " ms"
              << "   (control period dt = " << dt_ctrl * ms << " ms)\n"
              << "  worst per control step: " << t_solve_worst * ms << " ms\n"
              << "  function evaluations  : " << t_feval_total * ms << " ms ("
              << (t_solve_total > 0 ? 100. * t_feval_total / t_solve_total : 0.) << " %)\n"
              << "  search direction (KKT): " << t_sd_total * ms << " ms ("
              << (t_solve_total > 0 ? 100. * t_sd_total / t_solve_total : 0.) << " %)\n";
    const bool realtime = (t_solve_worst < dt_ctrl);
    std::cout << "  -> " << (realtime ? "real-time feasible: worst-case solve fits in dt."
                                      : "NOT real-time: worst-case solve exceeds dt.") << "\n";
    std::cout << "\nwrote quadrotor_mpc.csv (" << N_steps << " control steps)\n";

    const bool recovered = (final_p < 0.1 && final_tilt < 5.0 && final_om < 0.5);
    std::cout << (recovered ? "RECOVERED to upright hover.\n" : "did NOT fully settle.\n");
    return recovered ? 0 : 1;
}
