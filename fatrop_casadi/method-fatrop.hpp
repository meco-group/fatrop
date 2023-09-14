#pragma once
#include <casadi/casadi.hpp>
#include <vector>
#include <memory>
#include "ocp/StageOCPApplication.hpp"
#include "function-evaluation.hpp"
#include "stage-problem.hpp"
#include "fatrop-casadi-problem.hpp"
/**
 *  Wraps a SingleStageProblem into a FatropCasadiProblem
 */
namespace fatrop
{
    namespace fatrop_casadi
    {
        class StageProblemInternal;
        class StageProblemFatropMethod : public FatropCasadiProblem, public StageMethod
        {
        public:
            StageProblemFatropMethod(StageProblem *problem) : problem(problem)
            {
            }
            // transcribe should to the following tings
            // 1) define the micros stages and add them to the FatropCasadiProblem
            // 2) add gist to the StageProblem
            void transcribe(const int K)
            {
                K_ = K;
                x_vec = casadi::MX::veccat(problem->states);
                u_vec = casadi::MX::veccat(problem->controls);
                const int nx = x_vec.size1();
                const int nu = u_vec.size1();
                initial_syms = MicroStageSyms{casadi::MX::sym("x0", nx), casadi::MX::sym("u0", nu), casadi::MX::sym("p_stage", 0), casadi::MX::sym("p_global", 0)};
                middle_syms = MicroStageSyms{casadi::MX::sym("x", nx), casadi::MX::sym("u", nu), casadi::MX::sym("p_stage", 0), casadi::MX::sym("p_global", 0)};
                terminal_syms = MicroStageSyms{casadi::MX::sym("xf", nx), casadi::MX(), casadi::MX::sym("p_stage", 0), casadi::MX::sym("p_global", 0)};
                // initialize dynamics
                casadi::MX x_next_vec;
                for (auto &x : problem->states)
                {
                    x_next_vec = cs::MX::veccat({x_next_vec, problem->x_next[x]});
                }
                casadi::MX g_eq_initial;
                casadi::MX g_eq_middle;
                casadi::MX g_eq_terminal;
                casadi::MX g_ineq_initial;
                casadi::MX g_ineq_middle;
                casadi::MX g_ineq_terminal;
                casadi::DM lb_initial;
                casadi::DM lb_middle;
                casadi::DM lb_terminal;
                casadi::DM ub_initial;
                casadi::DM ub_middle;
                casadi::DM ub_terminal;
                // initialize constraints
                auto t_mode = MXPlaceholder::evaluation_mode::transcribe;
                for (auto &constraint : problem->constraints)
                {
                    auto placeholder = problem->ocp->placeholders(constraint, t_mode);
                    auto constraint_helped = ConstraintHelper(placeholder);
                    if (constraint.at_t0)
                    {

                        g_eq_initial = casadi::MX::veccat({g_eq_initial, problem->ocp->placeholders(problem->at_t0(constraint_helped.g), t_mode)});
                        g_ineq_initial = casadi::MX::veccat({g_ineq_initial, problem->ocp->placeholders(problem->at_t0(constraint_helped.g_ineq), t_mode)});
                        lb_initial = casadi::DM::veccat({lb_initial, constraint_helped.lb});
                        ub_initial = casadi::DM::veccat({ub_initial, constraint_helped.ub});
                    }
                    if (constraint.at_path)
                    {
                        g_eq_middle = casadi::MX::veccat({g_eq_middle, problem->ocp->placeholders(problem->at_path(constraint_helped.g), t_mode)});
                        g_ineq_middle = casadi::MX::veccat({g_ineq_middle, problem->ocp->placeholders(problem->at_path(constraint_helped.g_ineq), t_mode)});
                        lb_middle = casadi::DM::veccat({lb_middle, constraint_helped.lb});
                        ub_middle = casadi::DM::veccat({ub_middle, constraint_helped.ub});
                    }
                    if (constraint.at_tf)
                    {
                        // std::cout << "constraint" << std::endl;
                        // std::cout << constraint << std::endl;
                        // std::cout << "placeholder" << std::endl;
                        // std::cout << placeholder << std::endl;
                        // std::cout << "constraint_helped" << std::endl;
                        // std::cout << constraint_helped << std::endl;

                        g_eq_terminal = casadi::MX::veccat({g_eq_terminal, problem->ocp->placeholders(problem->at_tf(constraint_helped.g), t_mode)});
                        g_ineq_terminal = casadi::MX::veccat({g_ineq_terminal, problem->ocp->placeholders(problem->at_tf(constraint_helped.g_ineq), t_mode)});
                        lb_terminal = casadi::DM::veccat({lb_terminal, constraint_helped.lb});
                        ub_terminal = casadi::DM::veccat({ub_terminal, constraint_helped.ub});
                    }
                }
                // add the objective
                casadi::MX obj_t0 = 0;
                casadi::MX obj_path = 0;
                casadi::MX obj_tf = 0;
                // iterate over rpoblem->objective terms
                for (auto &objective : problem->objective_terms)
                {
                    // MXPlaceholder ph = problem->placeholders[objective];
                    if (objective.at_t0)
                    {
                        obj_t0 += problem->ocp->placeholders(objective, t_mode);
                    }
                    if (objective.at_path)
                    {
                        obj_path += problem->ocp->placeholders(objective, t_mode);
                    }
                    if (objective.at_tf)
                    {
                        obj_tf += problem->ocp->placeholders(objective, t_mode);
                    }
                }
                auto microstage_initial = MicroStage::create(initial_syms, obj_t0, problem->ocp->placeholders(problem->at_t0(x_next_vec), t_mode), g_eq_initial, g_ineq_initial, DM_to_vec_helper::DM_to_vec(lb_initial), DM_to_vec_helper::DM_to_vec(ub_initial));
                auto microstage_middle = MicroStage::create(middle_syms, obj_path, problem->ocp->placeholders(problem->at_path(x_next_vec), t_mode), g_eq_middle, g_ineq_middle, DM_to_vec_helper::DM_to_vec(lb_middle), DM_to_vec_helper::DM_to_vec(ub_middle));
                auto microstage_terminal = MicroStage::create(terminal_syms, obj_tf, cs::MX(), g_eq_terminal, g_ineq_terminal, DM_to_vec_helper::DM_to_vec(lb_terminal), DM_to_vec_helper::DM_to_vec(ub_terminal));

                // add the microstages
                push_back(microstage_initial);
                for (int i = 0; i < K - 1; i++)
                {
                    push_back(microstage_middle);
                }
                push_back(microstage_terminal);
            }

        protected:
            friend class Placeholders;
            casadi::MX eval_at(const PlaceHolderType type, const casadi::MX &expr)
            {
                switch (type)
                {
                case (PlaceHolderType::at_t0):
                    return casadi::MX::substitute({expr}, {x_vec, u_vec, p_stage_vec, p_global_vec}, {initial_syms.x, initial_syms.u, initial_syms.p_stage, initial_syms.p_global})[0];
                case (PlaceHolderType::at_tf):
                    return casadi::MX::substitute({expr}, {x_vec, u_vec, p_stage_vec, p_global_vec}, {terminal_syms.x, terminal_syms.u, terminal_syms.p_stage, terminal_syms.p_global})[0];
                case (PlaceHolderType::at_path):
                    return casadi::MX::substitute({expr}, {x_vec, u_vec, p_stage_vec, p_global_vec}, {middle_syms.x, middle_syms.u, middle_syms.p_stage, middle_syms.p_global})[0];
                default:
                    throw std::runtime_error("Unknown placeholder type");
                }
            }
            int K_;
            StageProblem *problem;
            casadi::MX x_vec;
            casadi::MX u_vec;
            casadi::MX p_stage_vec;
            casadi::MX p_global_vec;
            std::vector<casadi::MX> x_gist;
            std::vector<casadi::MX> u_gist;
            std::vector<casadi::MX> p_stage_gist;
            casadi::MX p_global_gist;
            MicroStageSyms initial_syms;
            MicroStageSyms middle_syms;
            MicroStageSyms terminal_syms;
        };
    }
};