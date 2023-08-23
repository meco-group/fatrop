#pragma once
#include <casadi/casadi.hpp>
#include "ocp/OCPAbstract.hpp"
#include <vector>
#include "shared-obj.hpp"
namespace fatrop
{
    namespace fatrop_casadi
    {
        namespace cs = casadi;
        struct MicroStageDims
        {
            int nx;
            int nxp1;
            int nu;
            int ng_equality;
            int ng_inequality;
            int np_stage;
            int np_global;
        };
        struct MicroStageSyms
        {
            cs::MX x;
            cs::MX u;
            cs::MX p_stage;
            cs::MX p_global;
        };
        struct MicroStageInternal
        {
            MicroStageInternal(const MicroStageDims &dims, const cs::Function &L, const cs::Function &dynamics, const cs::Function &g_eq, const cs::Function &g_ineq, const std::vector<double> &lb, const std::vector<double> &ub)
                : dims(dims), L(L), g_equality(g_eq), g_inequality(g_ineq), Lb(lb), Ub(ub)
            {
                // instantiate required symbols
                cs::MX x = cs::MX::sym("x", dims.nx);
                cs::MX xp1 = cs::MX::sym("xp1", dims.nxp1);
                cs::MX u = cs::MX::sym("u", dims.nu);
                cs::MX p_stage = cs::MX::sym("p", dims.np_stage);
                cs::MX p_global = cs::MX::sym("p_global", dims.np_global);
                // instantiate required dual variabels symbols
                cs::MX lam_g_equality = cs::MX::sym("lam_g_equality", dims.ng_equality);
                cs::MX lam_g_inequality = cs::MX::sym("lam_g_inequality", dims.ng_inequality);
                cs::MX lam_dynamics = cs::MX::sym("lam_dynamics", dims.nxp1);
                // instantiate all requried symbols
                auto L_sym = L({x, u, p_stage, p_global})[0];
                auto ux = cs::MX::vertcat({u, x});
                auto rq_sym = cs::MX();
                auto RSQ_sym = cs::MX::hessian(L_sym, ux, rq_sym);
                // deduce rq
                rq = casadi::Function("rq", {x, u, p_stage, p_global}, {cs::MX::densify(rq_sym)});
                // deduce BAbt
                if (dims.nxp1 > 0)
                {
                    auto x_next_sym = dynamics({x, u, p_stage, p_global})[0];
                    auto b_sym = -xp1 + x_next_sym;
                    b = casadi::Function("b", {x, xp1, u, p_stage, p_global}, {cs::MX::densify(b_sym)});
                    auto BAbt_sym = cs::MX::horzcat({cs::MX::jacobian(x_next_sym, ux), b_sym}).T();
                    BAbt = cs::Function("BAbt", {x, xp1, u, p_stage, p_global}, {cs::MX::densify(BAbt_sym)});
                    auto rq_dyn_sym = cs::MX();
                    auto RSQ_dyn_sym = cs::MX::hessian(cs::MX::dot(x_next_sym, lam_dynamics), ux, rq_dyn_sym);
                    rq_sym += rq_dyn_sym;
                    RSQ_sym += RSQ_dyn_sym;
                }
                // deduce Ggt_equality
                if (dims.ng_equality > 0)
                {
                    auto g_eq_sym = g_equality({x, u, p_stage, p_global})[0];
                    auto Ggt_equality_sym = cs::MX::horzcat({cs::MX::jacobian(g_eq_sym, ux), g_eq_sym}).T();
                    Ggt_equality = cs::Function("Ggt_equality", {x, u, p_stage, p_global}, {cs::MX::densify(Ggt_equality_sym)});
                    auto rq_eq = cs::MX();
                    auto RSQ_eq = cs::MX::hessian(cs::MX::dot(g_eq_sym, lam_g_equality), ux, rq_eq);
                    rq_sym += rq_eq;
                    RSQ_sym += RSQ_eq;
                }
                // deduce Ggt_inequality
                if (dims.ng_inequality > 0)
                {
                    auto g_ineq_sym = g_inequality({x, u, p_stage, p_global})[0];
                    auto Ggt_inequality_sym = cs::MX::horzcat({cs::MX::jacobian(g_ineq_sym, ux), g_ineq_sym}).T();
                    Ggt_inequality = cs::Function("Ggt_inequality", {x, u, p_stage, p_global}, {cs::MX::densify(Ggt_inequality_sym)});
                    auto rq_ineq = cs::MX();
                    auto RSQ_ineq = cs::MX::hessian(cs::MX::dot(g_ineq_sym, lam_g_inequality), ux, rq_ineq);
                    rq_sym += rq_ineq;
                    RSQ_sym += RSQ_ineq;
                }
                // deduce RSQrq
                RSQrq = casadi::Function("RSQrq", {x, u, p_stage, p_global, lam_dynamics, lam_g_equality, lam_g_inequality}, {cs::MX::densify(cs::MX::horzcat({RSQ_sym, rq_sym}).T())});
            }
            MicroStageInternal(const MicroStageSyms &syms, const cs::MX &L, const cs::MX &x_next, const cs::MX &g_eq, const cs::MX &g_ineq, const std::vector<double> &lb, const std::vector<double> &ub)
                : MicroStageInternal(MicroStageDims{(int)syms.x.size1(), (int)x_next.size1(), (int)syms.u.size1(), (int) g_eq.size1(), (int) g_ineq.size1(), (int)syms.p_stage.size1(), (int)syms.p_global.size1()}, cs::Function("L", {syms.x, syms.u, syms.p_stage, syms.p_global}, {L}), cs::Function("dynamics", {syms.x, syms.u, syms.p_stage, syms.p_global}, {x_next}), cs::Function("g_eq", {syms.x, syms.u, syms.p_stage, syms.p_global}, {g_eq}), cs::Function("g_ineq", {syms.x, syms.u, syms.p_stage, syms.p_global}, {g_ineq}), lb, ub)
            {
            }
            void expand()
            {
                // expand all functions
                L.expand();
                if (dims.nxp1 > 0)
                    b.expand();
                if (dims.ng_equality > 0)
                    g_equality.expand();
                if (dims.ng_inequality > 0)
                    g_inequality.expand();
                RSQrq.expand();
                rq.expand();
                if (dims.ng_equality > 0)
                    Ggt_equality.expand();
                if (dims.ng_inequality > 0)
                    Ggt_inequality.expand();
                if (dims.nxp1 > 0)
                    BAbt.expand();
            }
            // quantities that must be provided
            MicroStageDims dims;
            cs::Function L;
            cs::Function b; // actually provided as dynamics
            cs::Function g_equality;
            cs::Function g_inequality;
            std::vector<double> Lb;
            std::vector<double> Ub;

            // deducable quantities, todo: make it possible to provide these
            cs::Function RSQrq;
            cs::Function rq;
            cs::Function Ggt_equality;
            cs::Function Ggt_inequality;
            cs::Function BAbt;
        };

        struct MicroStage : public SharedObj<MicroStageInternal, MicroStage>
        {
            using SharedObj<MicroStageInternal, MicroStage>::SharedObj;
            MicroStage &expand()
            {
                (*this)->expand();
                return *this;
            }
        };

        class FatropCasadiProblem : public std::vector<MicroStage>
        {
        };
    } // namespace casadi
} // namespace fatrop
