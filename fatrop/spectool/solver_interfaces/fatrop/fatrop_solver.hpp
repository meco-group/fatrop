#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
#include "fatrop/spectool/solver_interfaces/solver.hpp"
#include "fatrop_function.hpp"
#include "fatrop_ocp_impl.hpp"
namespace fatrop
{
    namespace spectool
    {
        namespace cs = casadi;
        struct uStageEvaluator
        {
        };
        class SolverFatrop : public SolverInterface
        {
        public:
            void transcribe(const Ocp &ocp_, const cs::Dict &opts)
            {
                fatrop_impl = std::make_shared<FatropOcpImpl>(ocp_, opts);
            }
            cs::Function to_function(const std::string &name, const Ocp &ocp_, std::vector<cs::MX> &gist_in, std::vector<cs::MX> &gist_out, const cs::Dict &opts)
            {
                gist(ocp_, gist_in, gist_out);
                // return FatropFunction(name, fatrop_impl, opts);
                // C-api approach
                // C_api
                C_api_userdata* userdata = new C_api_userdata(fatrop_impl, opts); 
                // userdata->ref_count = 1;
                std::cout << "userdata pointer is " << userdata << std::endl;
                auto func = cs::external("fatrop_func", "/home/lander/fatrop/build/fatrop/spectool/libspectool.so", cs::Dict{{"user_data", static_cast<void*>(userdata)}});
                // fatrop_func_decref(static_cast<void*>(userdata));
                return func;

            };
            void gist(const Ocp &ocp_, std::vector<cs::MX> &in, std::vector<cs::MX> &out)
            {
                std::vector<cs::MX> variables_v;
                std::vector<cs::MX> control_grid_p_v;
                // add gist
                for (auto ustage_it = ocp_.get_ustages().begin(); ustage_it != ocp_.get_ustages().end(); ustage_it++)
                {
                    const auto &ustage = *ustage_it;
                    const auto &prev = ustage_it == ocp_.get_ustages().begin() ? std::make_shared<uStageInternal>() : (ustage_it - 1)->get_internal();
                    auto controls = cs::MX::veccat(ustage.get_controls(true, prev));
                    auto states = cs::MX::veccat(ustage.get_states(true, prev));
                    auto control_grid_p = cs::MX::veccat(ustage.get_control_parameters());
                    for (int k = 0; k < ustage.K(); k++)
                    {
                        variables_v.push_back(ustage.eval_at_control(controls, k));
                        variables_v.push_back(ustage.eval_at_control(states, k));
                        control_grid_p_v.push_back(ustage.eval_at_control(control_grid_p, k));
                    }
                }
                in = {cs::MX::veccat(variables_v), cs::MX::veccat(control_grid_p_v), cs::MX::veccat(ocp_.get_global_parameters())};
                out = {cs::MX::veccat(variables_v)};
            }
            std::shared_ptr<FatropOcpImpl> fatrop_impl;
        };
    } // namespace spectrop
} // namespace fatrop