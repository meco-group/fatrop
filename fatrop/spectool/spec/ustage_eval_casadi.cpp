#include "ustage_eval_casadi.hpp"

namespace fatrop
{
    namespace spectool
    {
        FatropuStageEvalCasadi::FatropuStageEvalCasadi(const uStageQuantities &sq, const cs::Dict &opts, CasadiJitCache &eval_cache)
        {
            auto optss = opts;
            bool expand = cs::get_from_dict(optss, "post_expand", true);
            bool jit = cs::get_from_dict(optss, "jit", false);
            // we handle jit ourself
            optss["jit"] = false;
            cs::Dict jit_options_ = cs::get_from_dict(optss, "jit_options", casadi::Dict({{"flags", "-Ofast -march=native -ffast-math"}}));
            K_ = sq.K;
            nu_ = nu(sq);
            nx_ = nx(sq);
            np_stage_ = np_stage(sq);
            np_global_ = np_global(sq);
            ng_eq_ = ng_eq(sq);
            ng_ineq_ = ng_ineq(sq);
            nxp1_ = nxp1(sq);
            cs::MX lam_g_equality = cs::MX::sym("lam_g_equality", ng_eq_);
            cs::MX lam_g_inequality = cs::MX::sym("lam_g_inequality", ng_ineq_);
            cs::MX lam_dynamics = cs::MX::sym("lam_dynamics", nxp1_);
            auto xp1 = cs::MX::sym("xp1", nxp1_);
            std::vector<cs::MX> ustage_syms{sq.u, sq.x, sq.p_stage, sq.p_global};

            auto hess_lag_obj = hess_lag_obj_sym(sq);
            auto dynamics_jac = dynamics_jacobian_sym(sq, xp1);
            auto equality_jac = equality_jacobian_sym(sq);
            auto inequality_jac = inequality_jacobian_sym(sq);
            auto hess_lag = hess_lag_sym(sq, lam_dynamics, lam_g_equality, lam_g_inequality);

            auto BAbt = cs::MX::horzcat({dynamics_jac.first, dynamics_jac.second}).T();
            auto Ggt_equality = cs::MX::horzcat({equality_jac.first, equality_jac.second}).T();
            auto Ggt_inequality = cs::MX::horzcat({inequality_jac.first, inequality_jac.second}).T();
            auto RSQrqt = cs::MX::horzcat({hess_lag.first, hess_lag.second}).T();

            L_ = CasadiFEWrap(cs::Function("L", ustage_syms, {sq.L}, optss), expand, jit, jit_options_, eval_cache);
            rq_ = CasadiFEWrap(cs::Function("rq", ustage_syms, {cs::MX::densify(hess_lag_obj.second)}, optss), expand, jit, jit_options_, eval_cache);
            b_ = CasadiFEWrap(cs::Function("b", {sq.x, xp1, sq.u, sq.p_stage, sq.p_global},
                                           {densify(dynamics_jac.second)}, optss),
                              expand, jit, jit_options_, eval_cache);
            BAbt_ = CasadiFEWrap(cs::Function("BAbt", {sq.x, xp1, sq.u, sq.p_stage, sq.p_global},
                                              {cs::MX::densify(BAbt)}, optss),
                                 expand, jit, jit_options_, eval_cache);
            Ggt_equality_ = CasadiFEWrap(cs::Function("Ggt_equality", ustage_syms, {cs::MX::densify(Ggt_equality)}, optss), expand, jit, jit_options_, eval_cache);
            g_equality_ = CasadiFEWrap(cs::Function("g_equality", ustage_syms, {cs::MX::densify(equality_jac.second)}, optss), expand, jit, jit_options_, eval_cache);
            Ggt_inequality_ = CasadiFEWrap(cs::Function("Ggt_inequality", ustage_syms, {cs::MX::densify(Ggt_inequality)}, optss), expand, jit, jit_options_, eval_cache);
            g_inequality_ = CasadiFEWrap(cs::Function("g_inequality", ustage_syms, {cs::MX::densify(inequality_jac.second)}, optss), expand, jit, jit_options_, eval_cache);
            RSQrqt_ = CasadiFEWrap(cs::Function("RSQrqt", {sq.u, sq.x, lam_dynamics, lam_g_equality, lam_g_inequality, sq.p_stage, sq.p_global}, {cs::MX::densify(RSQrqt)}, optss), expand, jit, jit_options_, eval_cache);
            Ub_ = sq.ub.nonzeros();
            Lb_ = sq.lb.nonzeros();
        }
        int FatropuStageEvalCasadi::nu(const uStageQuantities &sq)
        {
            return sq.nu();
        }
        int FatropuStageEvalCasadi::nx(const uStageQuantities &sq)
        {
            return sq.nx();
        }
        int FatropuStageEvalCasadi::np_stage(const uStageQuantities &sq)
        {
            return sq.np_stage();
        }
        int FatropuStageEvalCasadi::np_global(const uStageQuantities &sq)
        {
            return sq.np_global();
        }
        int FatropuStageEvalCasadi::ng_eq(const uStageQuantities &sq)
        {
            return sq.ng_eq();
        }
        int FatropuStageEvalCasadi::ng_ineq(const uStageQuantities &sq)
        {
            return sq.ng_ineq();
        }
        int FatropuStageEvalCasadi::nxp1(const uStageQuantities &sq)
        {
            return sq.nxp1();
        }
        std::pair<cs::MX, cs::MX> FatropuStageEvalCasadi::dynamics_jacobian_sym(const uStageQuantities &sq, const cs::MX &xp1)
        {
            auto ux = cs::MX::veccat({sq.u, sq.x});
            auto G = cs::MX::jacobian(sq.x_next, ux);
            auto g = -xp1 + sq.x_next;
            return std::make_pair(G, g);
        }

        std::pair<cs::MX, cs::MX> FatropuStageEvalCasadi::equality_jacobian_sym(const uStageQuantities &sq)
        {
            auto ux = cs::MX::veccat({sq.u, sq.x});
            auto G = cs::MX::jacobian(sq.g, ux);
            auto g = sq.g;
            return std::make_pair(G, g);
        }

        std::pair<cs::MX, cs::MX> FatropuStageEvalCasadi::inequality_jacobian_sym(const uStageQuantities &sq)
        {
            auto ux = cs::MX::veccat({sq.u, sq.x});
            auto G = cs::MX::jacobian(sq.g_ineq, ux);
            auto g = sq.g_ineq;
            return std::make_pair(G, g);
        }

        std::pair<cs::MX, cs::MX> FatropuStageEvalCasadi::hess_lag_obj_sym(const uStageQuantities &sq)
        {
            auto ux = cs::MX::veccat({sq.u, sq.x});
            auto rq_sym = cs::MX();
            auto RSQ_sym = cs::MX::hessian(sq.L, ux, rq_sym);
            return std::make_pair(RSQ_sym, rq_sym);
        }
        std::pair<cs::MX, cs::MX> FatropuStageEvalCasadi::hess_lag_dyn_sym(const uStageQuantities &sq, const cs::MX &lam_dyn)
        {
            auto ux = cs::MX::veccat({sq.u, sq.x});
            auto rq_sym = cs::MX();
            auto RSQ_sym = cs::MX::hessian(cs::MX::dot(sq.x_next, lam_dyn), ux, rq_sym);
            return std::make_pair(RSQ_sym, rq_sym);
        }
        std::pair<cs::MX, cs::MX> FatropuStageEvalCasadi::hess_lag_eq_sym(const uStageQuantities &sq, const cs::MX &lam_g_equality)
        {
            auto ux = cs::MX::veccat({sq.u, sq.x});
            auto rq_sym = cs::MX();
            auto RSQ_sym = cs::MX::hessian(cs::MX::dot(sq.g, lam_g_equality), ux, rq_sym);
            return std::make_pair(RSQ_sym, rq_sym);
        }
        std::pair<cs::MX, cs::MX> FatropuStageEvalCasadi::hess_lag_ineq_sym(const uStageQuantities &sq, const cs::MX &lam_g_inequality)
        {
            auto ux = cs::MX::veccat({sq.u, sq.x});
            auto rq_sym = cs::MX();
            auto RSQ_sym = cs::MX::hessian(cs::MX::dot(sq.g_ineq, lam_g_inequality), ux, rq_sym);
            return std::make_pair(RSQ_sym, rq_sym);
        }
        std::pair<cs::MX, cs::MX> FatropuStageEvalCasadi::hess_lag_sym(const uStageQuantities &sq, const cs::MX &lam_dyn, const cs::MX &lam_g_equality, const cs::MX &lam_g_inequality)
        {
            auto hess_obj = hess_lag_obj_sym(sq);
            auto hess_dyn = hess_lag_dyn_sym(sq, lam_dyn);
            auto hess_eq = hess_lag_eq_sym(sq, lam_g_equality);
            auto hess_ineq = hess_lag_ineq_sym(sq, lam_g_inequality);
            auto RSQ_sym = hess_obj.first + hess_dyn.first + hess_eq.first + hess_ineq.first;
            auto rq_sym = hess_obj.second + hess_dyn.second + hess_eq.second + hess_ineq.second;
            return std::make_pair(RSQ_sym, rq_sym);
        }

        int FatropuStageEvalCasadi::get_nx() const { return nx_; };
        int FatropuStageEvalCasadi::get_nu() const { return nu_; };
        int FatropuStageEvalCasadi::get_ng() const { return ng_eq_; };
        int FatropuStageEvalCasadi::get_n_stage_params() const { return np_stage_; };
        int FatropuStageEvalCasadi::get_ng_ineq() const { return ng_ineq_; };
        int FatropuStageEvalCasadi::eval_BAbt(
            const double *states_kp1,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            MAT *res)
        {
            const double *arg[5];
            arg[0] = states_k;
            arg[1] = states_kp1;
            arg[2] = inputs_k;
            arg[3] = stage_params_k;
            arg[4] = global_params;
            BAbt_.eval(arg, res);
            return 0;
        };
        int FatropuStageEvalCasadi::eval_RSQrqt(
            const double *objective_scale,
            const double *inputs_k,
            const double *states_k,
            const double *lam_dyn_k,
            const double *lam_eq_k,
            const double *lam_eq_ineq_k,
            const double *stage_params_k,
            const double *global_params,
            MAT *res)
        {
            const double *arg[7];
            arg[0] = inputs_k;
            arg[1] = states_k;
            arg[2] = lam_dyn_k;
            arg[3] = lam_eq_k;
            arg[4] = lam_eq_ineq_k;
            arg[5] = stage_params_k;
            arg[6] = global_params;
            RSQrqt_.eval(arg, res);
            return 0;
        };
        int FatropuStageEvalCasadi::eval_Ggt(
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            MAT *res)
        {
            const double *arg[4];
            arg[0] = inputs_k;
            arg[1] = states_k;
            arg[2] = stage_params_k;
            arg[3] = global_params;
            Ggt_equality_.eval(arg, res);
            return 0;
        };
        int FatropuStageEvalCasadi::eval_Ggt_ineq(
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            MAT *res)
        {
            const double *arg[4];
            arg[0] = inputs_k;
            arg[1] = states_k;
            arg[2] = stage_params_k;
            arg[3] = global_params;
            Ggt_inequality_.eval(arg, res);
            return 0;
        };
        int FatropuStageEvalCasadi::eval_b(
            const double *states_kp1,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            double *res)
        {
            const double *arg[5];
            arg[0] = states_k;
            arg[1] = states_kp1;
            arg[2] = inputs_k;
            arg[3] = stage_params_k;
            arg[4] = global_params;
            b_.eval(arg, res);
            return 0;
        }
        int FatropuStageEvalCasadi::eval_g(
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            double *res)
        {

            const double *arg[4];
            arg[0] = inputs_k;
            arg[1] = states_k;
            arg[2] = stage_params_k;
            arg[3] = global_params;
            g_equality_.eval(arg, res);
            return 0;
        }
        int FatropuStageEvalCasadi::eval_gineq(
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            double *res)
        {
            const double *arg[4];
            arg[0] = inputs_k;
            arg[1] = states_k;
            arg[2] = stage_params_k;
            arg[3] = global_params;
            g_inequality_.eval(arg, res);
            return 0;
        }
        int FatropuStageEvalCasadi::eval_rq(
            const double *objective_scale,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            double *res)
        {
            const double *arg[4];
            arg[0] = inputs_k;
            arg[1] = states_k;
            arg[2] = stage_params_k;
            arg[3] = global_params;
            rq_.eval(arg, res);
            return 0;
        }
        int FatropuStageEvalCasadi::eval_L(
            const double *objective_scale,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            double *res)
        {
            const double *arg[4];
            arg[0] = inputs_k;
            arg[1] = states_k;
            arg[2] = stage_params_k;
            arg[3] = global_params;
            L_.eval(arg, res);
            return 0;
        }
        int FatropuStageEvalCasadi::get_bounds(double *lower, double *upper) const
        {
            for (int i = 0; i < ng_ineq_; i++)
            {
                lower[i] = Lb_[i];
                upper[i] = Ub_[i];
            }
            return 0;
        }

    }
} // namespace fatrop