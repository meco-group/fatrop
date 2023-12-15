#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include <memory>
namespace fatrop
{
    namespace spectool
    {
        namespace cs = casadi;

        class uStageInternal; // forward declaration

        struct uStageQuantities
        {
            cs::MX x;
            cs::MX xp1; 
            cs::MX u;
            cs::MX lam_dyn; 
            cs::MX lam_g_ineq; 
            cs::MX lam_g_eq; 
            cs::MX p_stage;
            cs::MX p_global;
            cs::MX L;
            cs::DM lb;
            cs::DM ub;
            std::pair<cs::MX, cs::MX> Gg_dyn; 
            std::pair<cs::MX, cs::MX> Gg_ineq;
            std::pair<cs::MX, cs::MX> Gg_eq;
            std::pair<cs::MX, cs::MX> hess_obj;
            std::pair<cs::MX, cs::MX> hess_dyn;
            std::pair<cs::MX, cs::MX> hess_g_ineq;
            std::pair<cs::MX, cs::MX> hess_g_eq;
            int K;
            int nu() const { return u.size1(); };
            int nx() const { return x.size1(); };
            int nxp1() const { return xp1.size1(); };
            int np_stage() const { return p_stage.size1(); };
            int np_global() const { return p_global.size1(); };
            int ng_eq() const { return Gg_eq.first.size1(); };
            int ng_ineq() const { return Gg_ineq.first.size1(); };
            static std::pair<cs::MX, cs::MX> dynamics_jacobian_sym(const cs::MX& ux, const cs::MX &xp1, const cs::MX& x_next)
            {
                auto G = cs::MX::jacobian(x_next, ux);
                auto g = -xp1 + x_next;
                return std::make_pair(G, g);
            }

            static std::pair<cs::MX, cs::MX> equality_jacobian_sym(const cs::MX& ux, const cs::MX &g)
            {
                auto G = cs::MX::jacobian(g, ux);
                return std::make_pair(G, g);
            }

            static std::pair<cs::MX, cs::MX> inequality_jacobian_sym(const cs::MX& ux, const cs::MX &g_ineq)
            {
                auto G = cs::MX::jacobian(g_ineq, ux);
                return std::make_pair(G, g_ineq);
            }

            static std::pair<cs::MX, cs::MX> hess_lag_obj_sym(const cs::MX& ux, const cs::MX& L)
            {
                auto rq_sym = cs::MX();
                auto RSQ_sym = cs::MX::hessian(L, ux, rq_sym);
                return std::make_pair(RSQ_sym, rq_sym);
            }
            static std::pair<cs::MX, cs::MX> hess_lag_dyn_sym(const cs::MX& ux, const cs::MX &x_next, const cs::MX &lam_dyn)
            {
                auto rq_sym = cs::MX();
                auto RSQ_sym = cs::MX::hessian(cs::MX::dot(x_next, lam_dyn), ux, rq_sym);
                return std::make_pair(RSQ_sym, rq_sym);
            }
            static std::pair<cs::MX, cs::MX> hess_lag_eq_sym(const cs::MX& ux, const cs::MX& g, const cs::MX &lam_g_equality)
            {
                auto rq_sym = cs::MX();
                auto RSQ_sym = cs::MX::hessian(cs::MX::dot(g, lam_g_equality), ux, rq_sym);
                return std::make_pair(RSQ_sym, rq_sym);
            }
            static std::pair<cs::MX, cs::MX> hess_lag_ineq_sym(const cs::MX&ux, const cs::MX&g_ineq, const cs::MX &lam_g_inequality)
            {
                auto rq_sym = cs::MX();
                auto RSQ_sym = cs::MX::hessian(cs::MX::dot(g_ineq, lam_g_inequality), ux, rq_sym);
                return std::make_pair(RSQ_sym, rq_sym);
            }

            static uStageQuantities create(const std::shared_ptr<const uStageInternal> &ustage, const std::shared_ptr<const uStageInternal> &prev, const std::shared_ptr<const uStageInternal> &next, const std::vector<cs::MX> &global_parameter_syms);
        };
    } // namespace spectrop
} // namespace fatrop