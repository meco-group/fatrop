#pragma once
#include <casadi/casadi.hpp>
#include "ocp/OCPAbstract.hpp"
#include <vector>
namespace fatrop
{
    namespace specification
    {
        namespace cs = casadi;
        struct StageQuantities
        {
            StageQuantities(const cs::Function &L, const cs::Function &dynamics, const cs::Function &g_eq, const cs::Function &g_ineq, const std::vector<double> &lb, const std::vector<double> &ub)
                : L(L), g_equality(g_eq), g_inequality(g_ineq), Lb(lb), Ub(ub)
            {
                // TODO, make this part of the solver?
                // this is not required for opti for example
                // deduce problem dimensions
                nx = L.nnz_in(0);
                nu = L.nnz_in(1);
                np_stage = L.nnz_in(2);
                nxp1 = dynamics.nnz_in(1);
                ng_equality = g_equality.nnz_out(0);
                ng_inequality = g_inequality.nnz_out(0);
                np_global = L.nnz_in(3);
                // instantiate required symbols
                cs::MX x = cs::MX::sym("x", nx);
                cs::MX xp1 = cs::MX::sym("xp1", nxp1);
                cs::MX u = cs::MX::sym("u", nu);
                cs::MX p_stage = cs::MX::sym("p", np_stage);
                cs::MX p_global = cs::MX::sym("p_global", np_global);
                // instantiate required dual variabels symbols
                cs::MX lam_g_equality = cs::MX::sym("lam_g_equality", ng_equality);
                cs::MX lam_g_inequality = cs::MX::sym("lam_g_inequality", ng_inequality);
                cs::MX lam_dynamics = cs::MX::sym("lam_dynamics", nxp1);
                // instantiate all requried symbols
                auto x_next_sym = dynamics({x, u, p_stage, p_global})[0];
                auto g_eq_sym = g_equality({x, u, p_stage, p_global})[0];
                auto g_ineq_sym = g_inequality({x, u, p_stage, p_global})[0];
                auto L_sym = L({x, u, p_stage, p_global})[0];
                auto ux = cs::MX::vertcat({x, u});
                auto rq_sym = cs::MX();
                auto RSQ_sym = cs::MX::hessian(L_sym, ux, rq_sym);
                // deduce BAbt
                if (nxp1 > 0)
                {
                    auto b_sym = -xp1 + x_next_sym;
                    b = casadi::Function("b", {x, xp1, u, p_stage, p_global}, {b_sym});
                    auto BAbt_sym = cs::MX::horzcat({cs::MX::jacobian(x_next_sym, ux), b_sym}).T();
                    BAbt = cs::Function("BAbt", {x,xp1, u, xp1, p_stage, p_global}, {BAbt_sym});
                    auto rq_dyn_sym = cs::MX();
                    auto RSQ_dyn_sym = cs::MX::hessian(cs::MX::dot(x_next_sym, lam_dynamics), ux, rq_dyn_sym);
                    rq_sym += rq_dyn_sym;
                    RSQ_sym += RSQ_dyn_sym;
                }
                // deduce Ggt_equality
                if (ng_equality > 0)
                {
                    auto Ggt_equality_sym = cs::MX::horzcat({cs::MX::jacobian(g_eq_sym, ux), g_eq_sym}).T();
                    Ggt_equality = cs::Function("Ggt_equality", {x, u, p_stage, p_global}, {Ggt_equality_sym});
                    auto rq_eq = cs::MX();
                    auto RSQ_eq = cs::MX::hessian(cs::MX::dot(g_eq_sym, lam_g_equality), ux, rq_eq);
                    rq_sym += rq_eq;
                    RSQ_sym += RSQ_eq;
                }
                // deduce Ggt_inequality
                if (ng_inequality > 0)
                {
                    auto Ggt_inequality_sym = cs::MX::horzcat({cs::MX::jacobian(g_ineq_sym, ux), g_ineq_sym}).T();
                    Ggt_inequality = cs::Function("Ggt_inequality", {x, u, p_stage, p_global}, {Ggt_inequality_sym});
                    auto rq_ineq = cs::MX();
                    auto RSQ_ineq = cs::MX::hessian(cs::MX::dot(g_ineq_sym, lam_g_inequality), ux, rq_ineq);
                    rq_sym += rq_ineq;
                    RSQ_sym += RSQ_ineq;
                }
                // deduce rq
                rq = casadi::Function("rq", {x, u, p_stage, p_global, lam_dynamics, lam_g_equality, lam_g_inequality}, {rq_sym});
                // deduce RSQrq
                RSQrq = casadi::Function("RSQrq", {x, u, p_stage, p_global, lam_dynamics, lam_g_equality, lam_g_inequality}, {RSQ_sym});
            }
            void expand()
            {
                // expand all functions
                L.expand();
                b.expand();
                g_equality.expand();
                g_inequality.expand();
                RSQrq.expand();
                rq.expand();
                Ggt_equality.expand();
                Ggt_inequality.expand();
                BAbt.expand();
            }
            // quantities that must be provided
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
            // problem dimensions
            int ng_equality;
            int ng_inequality;
            int nx;
            int nxp1;
            int nu;
            int np_stage;
            int np_global;
        };

        class FatropCasadiProblem : public std::vector<std::unique_ptr<StageQuantities>>
        {
        };
    } // namespace specification
} // namespace fatrop
