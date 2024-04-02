#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
#include "fatrop/spectool/solver_interfaces/solver.hpp"
#include "fatrop_function.hpp"
// #include "fatrop_ocp_impl.hpp"
// #include "fatrop/ocp/CasadiCApiWrap.cpp.in"
#include "fatrop/spectool/spec/ustage_eval_casadi.hpp"
#include "fatrop/ocp/UStageOCPImpl.hpp"
#include "fatrop/spectool/spec/ocp.hpp"
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
            static std::vector<std::shared_ptr<FatropuStageEvalCasadi>> build(const Ocp &ocp, const cs::Dict &opts)
            {
                std::vector<std::shared_ptr<FatropuStageEvalCasadi>> ret;
                int horizon_length_ = 0;
                CasadiJitCache eval_cache;
                std::unordered_map<uStageInternal *, std::shared_ptr<FatropuStageEvalCasadi>> ustage_eval_cache;
                std::shared_ptr<FatropuStageEvalCasadi> evaluator = nullptr;
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
                // check that there are no constraints over multiple timesteps
                if (ocp_.get_multi_ustage_constraints().size() > 0)
                    throw std::runtime_error("Fatrop solver does not support constraints over different time steps.");
                int n_global_parameters_ = cs::MX::veccat(ocp_.get_global_parameters()).size1();
                fatrop_impl = std::make_shared<UStageOCPImpl<FatropuStageEvalCasadi>>(UStageEvalBuilder::build(ocp_, opts), n_global_parameters_);
                funct_opts = opts;
            }
            cs::Dict stats()
            {
                if (app)
                {
                    auto stats = app->get_stats();
                    cs::Dict ret;
                    ret["iter_count"] = stats.iterations_count;
                    ret["t_total"] = stats.time_total;
                    ret["t_FE"] = stats.eval_hess_time + stats.eval_jac_time + stats.eval_cv_time + stats.eval_grad_time + stats.eval_obj_time;
                    ret["t_init"] = stats.initialization_time;
                    ret["t_duinf"] = stats.duinf_time;
                    ret["t_compute_sd"] = stats.compute_sd_time;
                    ret["t_eval_hess"] = stats.eval_hess_time;
                    ret["t_eval_jac"] = stats.eval_jac_time;
                    ret["t_eval_cv"] = stats.eval_cv_time;
                    ret["t_eval_grad"] = stats.eval_grad_time;
                    ret["t_eval_obj"] = stats.eval_obj_time;
                    ret["success"] = stats.return_flag == 0;
                    return ret;
                }
                return {{}};
            }
            cs::Function to_function(const std::string &name, const Ocp &ocp_, std::vector<cs::MX> &gist_in, std::vector<cs::MX> &gist_out, const cs::Dict &opts)
            {
                gist(ocp_, gist_in, gist_out);
                // return FatropFunction(name, fatrop_impl, opts);
                // C-api approach
                // C_api
                app = std::make_shared<fatrop::OCPAbstractApplication>(fatrop_impl);
                for (auto opt : opts)
                {
                    if (opt.second.is_double())
                        app->set_option(opt.first, (double)opt.second);
                    else if (opt.second.is_int())
                        app->set_option(opt.first, (int)opt.second);
                    else if (opt.second.is_bool())
                        app->set_option(opt.first, (bool)opt.second);
                    else if (opt.second.is_string())
                        app->set_option(opt.first, (std::string)opt.second);
                }
                app->build();
                // go over the options and set
                return FatropFunction(name, app, opts, funct_opts);
            };
            std::shared_ptr<fatrop::OCPAbstractApplication> app;
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
            std::shared_ptr<UStageOCPImpl<FatropuStageEvalCasadi>> fatrop_impl;
            cs::Dict funct_opts;
        };
    } // namespace spectrop
} // namespace fatrop