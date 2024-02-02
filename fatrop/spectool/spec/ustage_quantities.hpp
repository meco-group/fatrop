#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include <memory>
#include <algorithm>
#include "custom_jacobian.hpp"
#include "custom_hessian.hpp"
#include "fatrop/spectool/auxiliary/constraint_helper.hpp"
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
            static std::pair<cs::MX, cs::MX> generate_jacobian(const cs::MX &x, const cs::MX &g, const Jacobian &jac);
            static std::pair<cs::MX, cs::MX> generate_hessian(const cs::MX &x, cs::MX &lam, const cs::MX &J, const cs::MX &g, const Hessian &hess);
            static std::pair<cs::MX, cs::MX> generate_jacobian(const cs::MX &x, const std::vector<cs::MX> &g, const uo_map_mx<Jacobian> &jac);
            static std::pair<cs::MX, cs::MX> generate_hessian(const cs::MX& x, cs::MX& lam, const cs::MX &J, const std::vector<cs::MX> &g, const uo_map_mx<Hessian> &hess);
            static uStageQuantities create(const std::shared_ptr<const uStageInternal> &ustage, const std::shared_ptr<const uStageInternal> &prev, const std::shared_ptr<const uStageInternal> &next, const std::vector<cs::MX> &global_parameter_syms);
        };
    } // namespace spectrop
} // namespace fatrop