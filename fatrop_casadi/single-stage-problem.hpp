#pragma once
#include <casadi/casadi.hpp>
#include <vector>
#include <map>
#include "fatrop-casadi-problem.hpp"
/**
 *  This file has a minimal casadi representation of a single stage problem.
 */
namespace fatrop
{
    namespace fatrop_casadi
    {
        struct StageProblemDimensions
        {
            int K;
            int nx;
            int nu;
            int np_stage;
            int np_global;
        };

        enum PlaceHolderType
        {
            at_t0,
            at_tf,
            sum_wfirst_wlast,
            sum_wofirst_wlast,
            sum_wfirst_woast,
            sum_wofirst_wolast,
        };
        class StageProblemInternal;
        struct MXPlaceholder : public casadi::MX
        {
            MXPlaceholder(const casadi::MX &expr, PlaceHolderType type, StageProblemInternal *stage) : casadi::MX(expr), type(type), stage(stage){};
            PlaceHolderType type;
            StageProblemInternal *stage; // TODO check if this can be done in a safer way
        };
        class StageProblemInternal;

        class Placeholders : public std::map<casadi::MX, MXPlaceholder>
        {
        public:
        private:
            std::vector<casadi::MX> get_all_placeholders(const casadi::MX &expr)
            {
                auto ret = std::vector<casadi::MX>();
                auto syms = casadi::MX::symvar(expr);
                for (auto &sym : syms)
                {
                    if (find(sym) != end())
                    {
                        ret.push_back(find(sym)->first);
                    }
                }
                return ret;
            }
        };

        struct StageProblem : public SharedObj<StageProblemInternal, StageProblem>
        {
            using SharedObj<StageProblemInternal, StageProblem>::SharedObj;
        };

        class StageProblemInternal
        {
        public:
            casadi::MX at_tf(const casadi::MX &expr)
            {
                return new_placeholder_expression(expr, PlaceHolderType::at_tf);
            }
            casadi::MX at_t0(const casadi::MX &expr)
            {
                return new_placeholder_expression(expr, PlaceHolderType::at_t0);
            }
            casadi::MX sum(const casadi::MX &expr, bool include_first = false, bool include_last = false)
            {
                if (include_first && include_last)
                    return new_placeholder_expression(expr, PlaceHolderType::sum_wfirst_wlast);
                if (include_first && !include_last)
                    return new_placeholder_expression(expr, PlaceHolderType::sum_wfirst_woast);
                if (!include_first && include_last)
                    return new_placeholder_expression(expr, PlaceHolderType::sum_wofirst_wlast);
                if (!include_first && !include_last)
                    return new_placeholder_expression(expr, PlaceHolderType::sum_wofirst_wolast);
            }
            void subject_to(const casadi::MX &expr, bool include_first = false, bool include_last = false)
            {
                // check if expr is a placeholder
                // if(placeholders.find() != placeholders.end())
                // {

                // }
                if (include_first)
                {
                    constraints_point.push_back(at_t0(expr));
                }
                if (include_last)
                {
                    constraints_point.push_back(at_tf(expr));
                }
                constraints.push_back(expr);
            }
            void add_objective(const casadi::MX &expr)
            {
                objective_terms.push_back(expr);
            }
            casadi::MX new_placeholder_expression(const casadi::MX &expr, PlaceHolderType type)
            {
                casadi::MX placeholder;
                placeholders[placeholder] = MXPlaceholder(expr, type, this);
                return placeholder;
            };
            StageProblem parent;
            StageProblem child;
            Placeholders placeholders;
            std::shared_ptr<FatropCasadiSolver> method;
            std::map<casadi::MX, casadi::MX> x_next;
            StageProblemDimensions dims;
            std::vector<casadi::MX> constraints;
            std::vector<casadi::MX> constraints_point;
            std::vector<casadi::MX> objective_terms;
            std::vector<casadi::MX> states;
            std::vector<casadi::MX> controls;
        };
    }
} // namespace fatrop_specification