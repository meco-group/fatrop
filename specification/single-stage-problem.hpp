#pragma once
#include <casadi/casadi.hpp>
#include <vector>
/**
 *  This file has a minimal casadi representation of a single stage problem.
 */
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
            std::vector<double> Lb;
            std::vector<double> Ub;
            int ng_equality;
            int ng_inequality;
        };
        struct SingleStageDimensions
        {
            int K;
            int nx;
            int nu;
            int np_stage;
            int np_global;
        };

        struct SingleStageProblem
        {
            SingleStageDimensions dims;
            StageQuantities initial;
            StageQuantities middle;
            StageQuantities terminal;
        };
    }
} // namespace fatrop_specification