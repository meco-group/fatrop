#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include <memory>
namespace fatrop
{
    namespace spectrop
    {
        namespace cs = casadi;

        class StageInternal; // forward declaration

        struct StageQuantities
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
            int nu;
            int nx;
            int nxp1;
            int np_stage;
            int np_global;
            int ng_eq;
            int ng_ineq;

            static StageQuantities create(const std::shared_ptr<StageInternal> &stage);
        };
    } // namespace spectrop
} // namespace fatrop