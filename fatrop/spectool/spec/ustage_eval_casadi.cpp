#include "ustage_eval_casadi.hpp"

namespace fatrop
{
    namespace spectool
    {
        FatropuStageEvalCasadi::FatropuStageEvalCasadi(const uStageQuantities &sq, const cs::Dict &opts, CasadiJitCache &eval_cache)
        {
            auto optss = opts;
            bool expand = cs::get_from_dict(optss, "expand", true);
            // delete the expand option from the dictionary
            optss.erase("expand");
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
            auto lam_dynamics = sq.lam_dyn;
            auto lam_g_equality = sq.lam_g_eq;
            auto lam_g_inequality = sq.lam_g_ineq;
            auto xp1 = sq.xp1;
            auto obj_scale = cs::MX::sym("obj_scale");
            std::vector<cs::MX> ustage_syms{sq.u, sq.x, sq.p_stage, sq.p_global};
            std::vector<cs::MX> os_ustage_syms{obj_scale, sq.u, sq.x, sq.p_stage, sq.p_global};

            auto dynamics_jac = sq.Gg_dyn;
            auto equality_jac = sq.Gg_eq;
            auto inequality_jac = sq.Gg_ineq;
            auto hess_lag_obj = sq.hess_obj;
            auto hess_lag_dyn = sq.hess_dyn;
            auto hess_lag_ineq = sq.hess_g_ineq;
            auto hess_lag_eq = sq.hess_g_eq;
            auto hess_lag = std::pair<cs::MX, cs::MX>{obj_scale*hess_lag_obj.first + hess_lag_dyn.first + hess_lag_ineq.first + hess_lag_eq.first,
                                                      obj_scale*hess_lag_obj.second + hess_lag_dyn.second + hess_lag_ineq.second + hess_lag_eq.second};

            auto BAbt = cs::MX::horzcat({dynamics_jac.first, dynamics_jac.second}).T();
            auto Ggt_equality = cs::MX::horzcat({equality_jac.first, equality_jac.second}).T();
            auto Ggt_inequality = cs::MX::horzcat({inequality_jac.first, inequality_jac.second}).T();
            auto RSQrqt = cs::MX::horzcat({hess_lag.first, hess_lag.second}).T();

            L_ = CasadiFEWrap(cs::Function("L", os_ustage_syms, {obj_scale*sq.L}, optss), expand, jit, jit_options_, eval_cache);
            rq_ = CasadiFEWrap(cs::Function("rq", os_ustage_syms, {cs::MX::densify(obj_scale*hess_lag_obj.second)}, optss), expand, jit, jit_options_, eval_cache);
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
            RSQrqt_ = CasadiFEWrap(cs::Function("RSQrqt", {obj_scale, sq.u, sq.x, lam_dynamics, lam_g_equality, lam_g_inequality, sq.p_stage, sq.p_global}, {cs::MX::densify(RSQrqt)}, optss), expand, jit, jit_options_, eval_cache);
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
            const double *arg[8];
            arg[0] = objective_scale;
            arg[1] = inputs_k;
            arg[2] = states_k;
            arg[3] = lam_dyn_k;
            arg[4] = lam_eq_k;
            arg[5] = lam_eq_ineq_k;
            arg[6] = stage_params_k;
            arg[7] = global_params;
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
            const double *arg[5];
            arg[0] = objective_scale;
            arg[1] = inputs_k;
            arg[2] = states_k;
            arg[3] = stage_params_k;
            arg[4] = global_params;
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
            const double *arg[5];
            arg[0] = objective_scale;
            arg[1] = inputs_k;
            arg[2] = states_k;
            arg[3] = stage_params_k;
            arg[4] = global_params;
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