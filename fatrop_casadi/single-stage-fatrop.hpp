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
        class StageProblemFatropMethod : public FatropCasadiProblem
        {
            StageProblemFatropMethod(const std::shared_ptr<StageProblemInternal> &problem) : problem(problem)
            {
            }
            void add_variables()
            {
            }
            void add_constraints()
            {
            }
            void add_objective()
            {
            }
            enum fill_placeholder_mode
            {
                transcribe,
                evaluate
            };

            casadi::MX fill_placeholder(const PlaceHolderType &type, const casadi::MX &expr, fill_placeholder_mode mode)
            {
                switch (mode)
                {
                case transcribe:
                    switch (type)
                    {
                    case at_t0:
                        return casadi::MX::substitute({expr}, {x_vec, u_vec, p_stage_vec, p_global_vec}, {initial_syms.x, initial_syms.u, initial_syms.p_stage, initial_syms.p_global})[0];
                    case at_tf:
                        return casadi::MX::substitute({expr}, {x_vec, u_vec, p_stage_vec, p_global_vec}, {terminal_syms.x, terminal_syms.u, terminal_syms.p_stage, terminal_syms.p_global})[0];
                    default:
                        throw std::runtime_error("Unknown placeholder type");
                    }
                case evaluate:
                    switch (type)
                    {
                    case at_t0:
                        return evaluate_at_control(expr, 0);
                    case at_tf:
                        return evaluate_at_control(expr, K - 1);
                    default:
                        throw std::runtime_error("Unknown placeholder type");
                    }
                }
            }
            casadi::MX evaluate_at_control(const casadi::MX &expr, const int k)
            {
                return expr;
            }
            int K;
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