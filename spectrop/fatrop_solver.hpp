#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "casadi_utilities.hpp"
#include "solver.hpp"
namespace fatrop
{
    namespace spectrop
    {
        namespace cs = casadi;
        struct StageEvaluator
        {

        };
        class SolverFatrop : public SolverInterface
        {
            void transcribe(const Ocp& ocp_) override;
        };
    } // namespace spectrop
} // namespace fatrop