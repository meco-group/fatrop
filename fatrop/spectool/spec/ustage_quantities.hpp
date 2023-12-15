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
            cs::MX u;
            cs::MX p_stage;
            cs::MX p_global;
            cs::MX L;
            cs::MX x_next;
            cs::MX g;
            cs::MX g_ineq;
            cs::DM lb;
            cs::DM ub;
            int K;
            int nu() const { return u.size1(); };
            int nx() const { return x.size1(); };
            int nxp1() const { return x_next.size1(); };
            int np_stage() const { return p_stage.size1(); };
            int np_global() const { return p_global.size1(); };
            int ng_eq() const { return g.size1(); };
            int ng_ineq() const { return g_ineq.size1(); };

            static uStageQuantities create(const std::shared_ptr<const uStageInternal> &ustage, const std::shared_ptr<const uStageInternal> &prev, const std::shared_ptr<const uStageInternal> &next, const std::vector<cs::MX> &global_parameter_syms);
        };
    } // namespace spectrop
} // namespace fatrop