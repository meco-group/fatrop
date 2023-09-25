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
                auto ret = std::make_shared<StageProblem>(*this);
                stageproblems.push_back(ret);
                return ret;
            }
            void subject_to(){

            };
            void solve()
            {
            }
            std::vector<std::shared_ptr<StageProblem>> stageproblems;
            std::vector<casadi::MX> x_gist;
        };

    } // namespace fatrop_casadi
} // namespace fatrop
