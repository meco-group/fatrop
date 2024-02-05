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
                auto app = std::make_shared<fatrop::OCPApplication>(fatrop_impl);
                app->build();
                // go over the options and set
                for (auto opt : opts)
                {
                    if (opt.second.is_double())
                        app->set_option(opt.first, (double)opt.second);
                    else if (opt.second.is_int())
                        app->set_option(opt.first, (int)opt.second);
                    else if (opt.second.is_bool())
                        app->set_option(opt.first, (bool)opt.second);
                }
                C_api_userdata *userdata = new C_api_userdata(app);
                userdata->ref_count = 0;
                // cs::Importer importer("/home/lander/fatrop/fatrop/ocp/liboldcapi.so", "dll");
                cs::Importer importer("/home/lander/fatrop/fatrop/ocp/Casadi_CApiWrap.cpp", "shell");
                reinterpret_cast<void (*)(C_api_userdata*)>(importer.get_function("set_user_data"))(userdata);
                auto func = cs::external("casadi_old_capi", importer);
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