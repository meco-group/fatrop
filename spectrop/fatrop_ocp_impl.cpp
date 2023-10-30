#include "fatrop_ocp_impl.hpp"
namespace fatrop{
    namespace spectrop
    {
            FatropStageEval FatropStageEval::create(const Stage &stage)
            {
                auto sq = StageQuantities::create(stage.get_internal());
                auto ret = FatropStageEval();
                ret.K = sq.K;
                ret.nu = sq.nu;
                ret.nx = sq.nx;
                ret.np_stage = sq.np_stage;
                ret.np_global = sq.np_global;
                ret.ng_eq = sq.ng_eq;
                ret.ng_ineq = sq.ng_ineq;
                cs::MX lam_g_equality = cs::MX::sym("lam_g_equality", sq.ng_eq);
                cs::MX lam_g_inequality = cs::MX::sym("lam_g_inequality", sq.ng_ineq);
                cs::MX lam_dynamics = cs::MX::sym("lam_dynamics", sq.nxp1);
                std::vector<cs::MX> stage_syms{sq.u, sq.x, sq.p_stage, sq.p_global};
                auto ux = cs::MX::veccat({sq.u, sq.x});
                ret.L = CasadiFEWrap(cs::Function("L", stage_syms, {sq.L}));
                auto rq_sym = cs::MX();
                auto RSQ_sym = cs::MX::hessian(sq.L, ux, rq_sym);
                ret.rq = CasadiFEWrap(cs::Function("rq", stage_syms, {cs::MX::densify(rq_sym)}));
                // if (sq.nxp1 > 0)
                {
                    auto xp1 = cs::MX::sym("xp1", sq.nxp1);
                    auto b = -xp1 + sq.x_next;
                    auto BAbt = cs::MX::horzcat({cs::MX::jacobian(sq.x_next, ux), b}).T();
                    ret.b = CasadiFEWrap(cs::Function("b", {sq.x, xp1, sq.u, sq.p_stage, sq.p_global},
                                                      {densify(b)}));
                    ret.BAbt = CasadiFEWrap(cs::Function("BAbt", {sq.x, xp1, sq.u, sq.p_stage, sq.p_global},
                                                         {cs::MX::densify(BAbt)}));
                    auto rq_dyn_sym = cs::MX();
                    auto RSQ_dyn_sym = cs::MX::hessian(cs::MX::dot(sq.x_next, lam_dynamics), ux, rq_dyn_sym);
                    rq_sym += rq_dyn_sym;
                    RSQ_sym += RSQ_dyn_sym;
                }
                // if(sq.ng_eq > 0)
                {
                    auto g_eq = sq.g;
                    auto Ggt_equality = cs::MX::horzcat({cs::MX::jacobian(g_eq, ux), g_eq}).T();
                    ret.Ggt_equality = CasadiFEWrap(cs::Function("Ggt_equality", stage_syms, {cs::MX::densify(Ggt_equality)}));
                    ret.g_equality = CasadiFEWrap(cs::Function("g_equality", stage_syms, {cs::MX::densify(g_eq)}));
                    auto rq_eq = cs::MX();
                    auto RSQ_eq = cs::MX::hessian(cs::MX::dot(g_eq, lam_g_equality), ux, rq_eq);
                    rq_sym += rq_eq;
                    RSQ_sym += RSQ_eq;
                }
                // if(sq.ng_ineq>0)
                {
                    auto g_ineq = sq.g_ineq;
                    auto Ggt_inequality = cs::MX::horzcat({cs::MX::jacobian(g_ineq, ux), g_ineq}).T();
                    ret.Ggt_inequality = CasadiFEWrap(cs::Function("Ggt_inequality", stage_syms, {cs::MX::densify(Ggt_inequality)}));
                    ret.g_inequality = CasadiFEWrap(cs::Function("g_inequality", stage_syms, {cs::MX::densify(g_ineq)}));
                    auto rq_ineq = cs::MX();
                    auto RSQ_ineq = cs::MX::hessian(cs::MX::dot(g_ineq, lam_g_inequality), ux, rq_ineq);
                    rq_sym += rq_ineq;
                    RSQ_sym += RSQ_ineq;
                }
                ret.RSQrqt = CasadiFEWrap(cs::Function("RSQrqt", {sq.u, sq.x, lam_dynamics, lam_g_equality, lam_g_inequality, sq.p_stage, sq.p_global}, {cs::MX::densify(cs::MX::horzcat({RSQ_sym, rq_sym}).T())}));
                ret.Ub = sq.ub.nonzeros();
                ret.Lb = sq.lb.nonzeros();
                return ret;
            }
    }
}