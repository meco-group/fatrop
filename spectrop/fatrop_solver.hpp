#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "casadi_utilities.hpp"
#include "solver.hpp"
#include "fatrop_function.hpp"
#include "fatrop_ocp_impl.hpp"
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
        public:
            void transcribe(const Ocp &ocp_) 
            {
                fatrop_impl = std::make_shared<FatropOcpImpl>(ocp_);
            }
            cs::Function to_function(const Ocp &ocp_, std::vector<cs::MX> &gist_in, std::vector<cs::MX> &gist_out)
            {
                std::vector<cs::MX> variables_v;
                std::vector<cs::MX> control_grid_p_v;
                // add gist 
                for(const auto& stage: ocp_.get_stages())
                {
                   variables_v.push_back(stage.sample(cs::MX::veccat(stage.get_controls())));
                   variables_v.push_back(stage.sample(cs::MX::veccat(stage.get_states())));
                   control_grid_p_v.push_back(stage.sample(cs::MX::veccat(stage.get_control_parameters())));
                }
                gist_in = {cs::MX::veccat(variables_v), cs::MX::veccat(control_grid_p_v), cs::MX::veccat(ocp_.get_global_parameters())};
                gist_out = {cs::MX::veccat(variables_v)};

                return FatropFunction(fatrop_impl);
            };
            std::shared_ptr<FatropOcpImpl> fatrop_impl;
        };
    } // namespace spectrop
} // namespace fatrop