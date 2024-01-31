#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
#include "fatrop/spectool/solver_interfaces/solver.hpp"
#include "fatrop/spectool/function_evaluation/casadi_fe.hpp"
#include "fatrop/spectool/spec/ustage.hpp"
#include "fatrop/spectool/spec/ocp.hpp"
#include "fatrop/spectool/spec/ustage_quantities.hpp"
#include "fatrop_function.hpp"
#include "fatrop/ocp/UStageOCPImpl.hpp"
namespace fatrop
{
    namespace spectool
    {
        namespace cs = casadi;
        struct uStageEvaluator
        {
        };
        struct UStageEvalBuilder
        {
            static std::vector<std::shared_ptr<UStageEvalAbstract>> build(const Ocp &ocp, const cs::Dict &opts)
            {
                std::vector<std::shared_ptr<UStageEvalAbstract>> ret; 
                int horizon_length_ = 0;
                CasadiJitCache eval_cache;
                std::unordered_map<uStageInternal *, std::shared_ptr<UStageEvalAbstract>> ustage_eval_cache;
                std::shared_ptr<UStageEvalAbstract> evaluator = nullptr;
                for (auto ustage_it = ocp.get_ustages().begin(); ustage_it != ocp.get_ustages().end(); ustage_it++)
                {
                    const auto &ustage = *ustage_it;
                    auto &original = ustage.get_original();
                    // check if evaluator is already available
                    auto ustage_eval_it = ustage_eval_cache.find(original.get());
                    if (original && ustage_eval_it != ustage_eval_cache.end())
                        evaluator = ustage_eval_it->second;
                    else
                    {
                        const auto &prev = ustage_it == ocp.get_ustages().begin() ? std::make_shared<uStageInternal>() : (ustage_it - 1)->get_internal();
                        const auto &next = ustage_it + 1 == ocp.get_ustages().end() ? nullptr : (ustage_it + 1)->get_internal();
                        evaluator = ustage.get_evaluator(prev, next, ocp.get_global_parameters(), opts, eval_cache);
                        // add to cache
                        ustage_eval_cache[ustage_it->get_original().get()] = evaluator;
                    }
                    ret.push_back(evaluator);
                    for (int i = 1; i < ustage.K(); i++)
                        ret.push_back(ret.back());
                    horizon_length_ += ustage.K();
                }
                return ret;
            }
        };


        class SolverFatrop : public SolverInterface
        {
        public:
            void transcribe(const Ocp &ocp_, const cs::Dict &opts)
            {
                int n_global_parameters_ = cs::MX::veccat(ocp_.get_global_parameters()).size1();
                fatrop_impl = std::make_shared<UStageOCPImpl>(UStageEvalBuilder::build(ocp_, opts), n_global_parameters_);
            }
            cs::Function to_function(const std::string &name, const Ocp &ocp_, std::vector<cs::MX> &gist_in, std::vector<cs::MX> &gist_out, const cs::Dict &opts)
            {
                gist(ocp_, gist_in, gist_out);
                return FatropFunction(name, fatrop_impl, opts);
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
            std::shared_ptr<UStageOCPImpl> fatrop_impl;
        };
    } // namespace spectrop
} // namespace fatrop