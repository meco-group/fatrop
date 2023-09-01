#pragma once
#include <casadi/casadi.hpp>
#include <vector>
#include <memory>
#include "ocp/StageOCPApplication.hpp"
#include "function-evaluation.hpp"
#include "single-stage-problem.hpp"
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
            StageProblemFatropMethod(const std::shared_ptr<StageProblemInternal> &problem) : problem(problem)
            {
            }
            void transcribe(const int K)
            {
                K_ = K;
                x_vec = casadi::MX::veccat(problem->states);
                u_vec = casadi::MX::veccat(problem->controls);
                const int nx = x_vec.size1();
                const int nu = u_vec.size1();
                initial_syms = MicroStageSyms{casadi::MX::sym("x0", nx), casadi::MX::sym("u0", nu), casadi::MX::sym("p_stage", 0), casadi::MX::sym("p_global", 0)};
                middle_syms = MicroStageSyms{casadi::MX::sym("x", nx), casadi::MX::sym("u", nu), casadi::MX::sym("p_stage", 0), casadi::MX::sym("p_global", 0)};
                terminal_syms = MicroStageSyms{casadi::MX::sym("xf", nx), casadi::MX::sym("uf", nu), casadi::MX::sym("p_stage", 0), casadi::MX::sym("p_global", 0)};
                // initialize dynamics
                casadi::MX x_next_vec;
                for (auto &x : problem->states)
                {
                    // assert(problem->x_next[x]);
                    x_next_vec.veccat({x_next_vec, problem->x_next[x]});
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
                for (auto &constraint : problem->constraints)
                {
                    bool at_t0;
                    bool at_tf;
                    bool at_path;
                    MXPlaceholder placeholder = problem->placeholders[constraint];
                    // assert(placeholder);
                    switch (placeholder.type)
                    {
                    case PlaceHolderType::at_t0:
                        at_t0 = true;
                        at_path = false;
                        at_tf = false;
                        break;
                    case PlaceHolderType::at_tf:
                        at_t0 = false;
                        at_path = false;
                        at_tf = true;
                        break;
                    case PlaceHolderType::path:
                        at_t0 = false;
                        at_path = true;
                        at_tf = false;
                        break;
                    case PlaceHolderType::path_t0tf:
                        at_t0 = true;
                        at_path = true;
                        at_tf = true;
                        break;
                    case PlaceHolderType::path_t0:
                        at_t0 = true;
                        at_path = true;
                        at_tf = false;
                        break;
                    case PlaceHolderType::path_tf:
                        at_t0 = false;
                        at_path = true;
                        at_tf = true;
                        break;
                    default:
                        throw std::runtime_error("Unknown placeholder type");
                    }
                    auto canon = casadi::Opti().advanced().canon_expr(constraint);
                }
            }

        private:
            void add_variables()
            {
            }
            void add_constraints()
            {
            }
            void add_objective()
            {
            }

        protected:
            friend class Placeholders;
            casadi::MX fill_placeholder(const PlaceHolderType &type, const casadi::MX &expr, MXPlaceholder::evaluation_mode mode)
            {

                switch (type)
                {
                case at_t0:
                    return evaluate_at_control(expr, 0, mode);
                case at_tf:
                    return evaluate_at_control(expr, K_ - 1, mode);
                case (path): //  && mode == MXPlaceholder::evaluation_mode::transcribe
                    return evaluate_at_control(expr, 1, mode);
                default:
                    throw std::runtime_error("Unknown placeholder type");
                }
            }

        private:
            casadi::MX evaluate_at_control(const casadi::MX &expr, const int k, MXPlaceholder::evaluation_mode mode)
            {
                switch (mode)
                {
                case MXPlaceholder::evaluation_mode::transcribe:
                {
                    if (k == 0)
                        return casadi::MX::substitute({expr}, {x_vec, u_vec, p_stage_vec, p_global_vec}, {initial_syms.x, initial_syms.u, initial_syms.p_stage, initial_syms.p_global})[0];
                    else if (k == K_ - 1)
                        return casadi::MX::substitute({expr}, {x_vec, u_vec, p_stage_vec, p_global_vec}, {terminal_syms.x, terminal_syms.u, terminal_syms.p_stage, terminal_syms.p_global})[0];
                    else
                        return casadi::MX::substitute({expr}, {x_vec, u_vec, p_stage_vec, p_global_vec}, {middle_syms.x, middle_syms.u, middle_syms.p_stage, middle_syms.p_global})[0];
                }
                case MXPlaceholder::evaluation_mode::evaluate:
                {
                    return casadi::MX::substitute({expr}, {x_vec, u_vec, p_stage_vec, p_global_vec}, {x_gist[k], u_gist[k], p_stage_gist[k], p_global_gist})[0];
                }
                default:
                {
                    throw std::runtime_error("Unknown evaluation type");
                }
                }
            }
            int K_;
            std::shared_ptr<StageProblemInternal> problem;
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
            MicroStage microstage_initial;
            MicroStage microstage_middle;
            MicroStage microstage_terminal;
        };
    }
};