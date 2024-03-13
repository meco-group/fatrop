#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
#include "fatrop/spectool/solver_interfaces/solver.hpp"
#include "fatrop/spectool/spec/ocp.hpp"
#include "fatrop/spectool/spec/ustage_quantities.hpp"

namespace fatrop
{
    namespace spectool
    {
        namespace cs = casadi;
        class SolverOpti : public SolverInterface
        {
        public:
            void transcribe(const Ocp &ocp_, const cs::Dict &opts)
            {
                optss_ = opts;
            }
            cs::Dict stats() 
            {
                return opti_function_.stats();
            }
            cs::Function to_function(const std::string &name, const Ocp &ocp, std::vector<cs::MX> &gist_in, std::vector<cs::MX> &gist_out, const cs::Dict &opts)
            {
                gist(ocp, gist_in, gist_out);
                int horizon_length_ = 0;
                int n_global_parameters_ = cs::MX::veccat(ocp.get_global_parameters()).size1();
                std::vector<uStageQuantities> ustages_quantities_;
                for (auto ustage_it = ocp.get_ustages().begin(); ustage_it != ocp.get_ustages().end(); ustage_it++)
                {
                    const auto &ustage = *ustage_it;
                    // auto &original = ustage.get_original();
                    const auto &prev = ustage_it == ocp.get_ustages().begin() ? std::make_shared<uStageInternal>() : (ustage_it - 1)->get_internal();
                    const auto &next = ustage_it + 1 == ocp.get_ustages().end() ? nullptr : (ustage_it + 1)->get_internal();
                    ustages_quantities_.push_back(uStageQuantities::create(ustage.get_internal(), prev, next, ocp.get_global_parameters()));
                    for (int i = 1; i < ustage.K(); i++)
                        ustages_quantities_.push_back(ustages_quantities_.back());
                    horizon_length_ += ustage.K();
                }
                cs::Opti opti;
                std::vector<cs::MX> variables_u;
                std::vector<cs::MX> variables_x;
                std::vector<cs::MX> variables_ux;
                std::vector<cs::MX> p_stage;
                // add all variables
                auto p_global = opti.parameter(n_global_parameters_);
                for (auto &ustage_q : ustages_quantities_)
                {
                    variables_u.push_back(opti.variable(ustage_q.nu()));
                    variables_x.push_back(opti.variable(ustage_q.nx()));
                    variables_ux.push_back(cs::MX::veccat({variables_u.back(), variables_x.back()}));
                    p_stage.push_back(opti.parameter(ustage_q.np_stage()));
                }
                // add dynamics constraints
                for (unsigned long i = 0; i < ustages_quantities_.size() - 1; i++)
                {
                    std::vector<cs::MX> from{ustages_quantities_[i].x, ustages_quantities_[i].u, ustages_quantities_[i].p_stage, ustages_quantities_[i].p_global, ustages_quantities_[i].xp1};
                    std::vector<cs::MX> to{variables_x[i], variables_u[i], p_stage[i], p_global, variables_x[i + 1]};
                    if (ustages_quantities_[i].nxp1() > 0)
                        opti.subject_to(cs::MX::substitute({ustages_quantities_[i].Gg_dyn.second == 0}, from, to)[0]);
                }
                // add path constraints
                for (unsigned long i = 0; i < ustages_quantities_.size(); i++)
                {
                    std::vector<cs::MX> from{ustages_quantities_[i].x, ustages_quantities_[i].u, ustages_quantities_[i].p_stage, ustages_quantities_[i].p_global};
                    std::vector<cs::MX> to{variables_x[i], variables_u[i], p_stage[i], p_global};
                    if (ustages_quantities_[i].ng_ineq() > 0)
                        opti.subject_to(cs::MX::substitute({ustages_quantities_[i].lb < ustages_quantities_[i].Gg_ineq.second <= ustages_quantities_[i].ub}, from, to)[0]);
                    if (ustages_quantities_[i].ng_eq())
                        opti.subject_to(cs::MX::substitute({ustages_quantities_[i].Gg_eq.second == 0}, from, to)[0]);
                }
                auto flat_variables_ux = cs::MX::veccat(variables_ux);
                auto flat_gist_out = cs::MX::veccat(gist_out);
                // add the constraints over multiple time steps
                for(auto& constr: ocp.get_multi_ustage_constraints())
                    opti.subject_to(cs::MX::substitute(constr, flat_gist_out,flat_variables_ux));
                // add objective
                cs::MX J = 0;
                for (unsigned long i = 0; i < ustages_quantities_.size(); i++)
                {
                    std::vector<cs::MX> from{ustages_quantities_[i].x, ustages_quantities_[i].u, ustages_quantities_[i].p_stage, ustages_quantities_[i].p_global};
                    std::vector<cs::MX> to{variables_x[i], variables_u[i], p_stage[i], p_global};
                    J += cs::MX::substitute({ustages_quantities_[i].L}, from, to)[0];
                }
                opti.minimize(J);
                opti.solver("ipopt", optss_, opts);
                auto all_vars = cs::MX::veccat({variables_ux});
                opti_function_ = opti.to_function(name, {all_vars, cs::MX::veccat(p_stage), p_global}, {all_vars});
                return opti_function_;
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
            cs::Dict optss_;
            cs::Function opti_function_;
        };
    } // namespace spectrop
} // namespace fatrop