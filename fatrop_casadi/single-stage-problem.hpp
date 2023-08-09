#pragma once
#include <casadi/casadi.hpp>
#include <vector>
#include "fatrop-casadi-problem.hpp"
/**
 *  This file has a minimal casadi representation of a single stage problem.
 */
namespace fatrop
{
    namespace fatrop_casadi
    {
        struct SingleStageDimensions
        {
            int K;
            int nx;
            int nu;
        };

        struct SingleStageProblem
        {
            SingleStageDimensions dims;
            MicroStageInternal initial;
            MicroStageInternal middle;
            MicroStageInternal terminal;
        };
    }
} // namespace fatrop_specification