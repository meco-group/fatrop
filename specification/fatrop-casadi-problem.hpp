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
            cs::Function L;
            cs::Function RSQrq;
            cs::Function rq;
            cs::Function g_equality;
            cs::Function g_inequality;
            cs::Function Ggt_equality;
            cs::Function Ggt_inequality;
            cs::Function BAbt;
            cs::Function b;
            std::vector<double> Lb;
            std::vector<double> Ub;
            int ng_equality;
            int ng_inequality;
            int nx;
            int nu;
            int np_stage;
        };

        class FatropCasadiProblem : public std::vector<StageQuantities>
        {
        };
    } // namespace specification
} // namespace fatrop
