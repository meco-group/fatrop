#pragma once
#include "stage-problem.hpp"
#include <vector>
#include <memory>
#include "ocp-internal.hpp"
namespace fatrop
{
    namespace fatrop_casadi
    {
        class Ocp : public std::shared_ptr<OcpInternal>
        {
        public:
            Ocp() : std::shared_ptr<OcpInternal>(std::make_shared<OcpInternal>())
            {
            }
            std::shared_ptr<StageProblem> stage()
            {
                return std::make_shared<StageProblem>(*this);
            }
            void subject_to(){

            };
            std::vector<std::shared_ptr<StageProblem>> stageproblems;
        };

    } // namespace fatrop_casadi
} // namespace fatrop
