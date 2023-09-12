#pragma once
#include "stage-problem.hpp"
#include <vector>
#include <memory>
namespace fatrop
{
    namespace fatrop_casadi
    {
        class Ocp
        {
            std::shared_ptr<StageProblem> stage()
            {
            }
            void subject_to()
            {

            };
            std::vector<std::shared_ptr<StageProblem>> stageproblems;
        };

    } // namespace fatrop_casadi
} // namespace fatrop
