#pragma once
#include <casadi/casadi.hpp>
#include <vector>
// #include <utility>
#include <functional>
#include <string>
#include <map>
#include "fatrop-casadi-problem.hpp"
#include "utilities.hpp"
#include "placeholders.hpp"
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
        
        enum class PlaceHolderType;

        class StageMethod
        {
        public:
            virtual casadi::MX fill_placeholder(const PlaceHolderType type, const casadi::MX &expr, MXPlaceholder::evaluation_mode mode)
            {
                throw std::runtime_error("No method set");
                return casadi::MX();
            };
            virtual void transcribe(const int K)
            {
                std::runtime_error("No method set");
            };
        };

        // struct StageProblem : public SharedObj<StageProblemInternal, StageProblem>
        // {
        //     using SharedObj<StageProblemInternal, StageProblem>::SharedObj;
        // };

        class StageProblem
        {
        public:
            casadi::MX state(const int m = 1, const int n = 1)
            {
                casadi::MX ret = casadi::MX::sym("x", m, n);
                states.push_back(ret);
                // x_next[ret] = casadi::MX::sym("dummy", 0, 0);
                return ret;
            }
            casadi::MX control(const int m = 1, const int n = 1)
            {
                casadi::MX ret = casadi::MX::sym("u", m, n);
                controls.push_back(ret);
                return ret;
            }
            void set_next(const casadi::MX &x, const casadi::MX &x_next_in)
            {
                this->x_next[x] = x_next_in;
            }
            casadi::MX at_tf(const casadi::MX &expr)
            {
                return new_placeholder_expression(expr, PlaceHolderType::at_tf, "at_tf");
            }
            casadi::MX at_t0(const casadi::MX &expr)
            {
                return new_placeholder_expression(expr, PlaceHolderType::at_t0, "at_t0");
            }
            casadi::MX at_path(const casadi::MX &expr)
            {
                return new_placeholder_expression(expr, PlaceHolderType::at_path, "at_path");
            }
            casadi::MX sum(const casadi::MX &expr, bool include_first = false, bool include_last = false)
            {
                casadi::MX sum = at_path(expr);
                if (include_first)
                    sum += at_t0(expr);
                if (include_last)
                    sum += at_tf(expr);
                return sum;
            }
            void subject_to(const casadi::MX &expr, bool include_first = false, bool include_last = false)
            {
                if (placeholders.has_placeholders(expr, {PlaceHolderType::at_t0})) // check for point constraints
                {
                    constraints.push_back(StageMX(expr, true, false, false));
                    return;
                }
                if (placeholders.has_placeholders(expr, {PlaceHolderType::at_tf})) // check for point constraints
                {
                    constraints.push_back(StageMX(expr, false, false, true));
                    return;
                }
                // path constraints
                constraints.push_back(StageMX(expr, include_first, true, include_last));
            }
            void add_objective(const casadi::MX &expr)
            {
                // get all terms of expr
                std::vector<casadi::MX> terms = get_terms_helper::get_terms(expr);
                // iterate over all terms
                for (auto &term : terms)
                {
                    // check the placeholder type
                    auto ph = placeholders.get_all_placeholders(term);
                    // check if only one placeholder
                    if (ph.size() != 1)
                        throw std::runtime_error("Objective term must contain exactly one placeholder");
                    auto type = ph[0].first.type;
                    objective_terms.push_back(StageMX(term, type == PlaceHolderType::at_t0, type == PlaceHolderType::at_path, type == PlaceHolderType::at_tf));
                }
            }
            casadi::MX fill_placeholder(const PlaceHolderType type, const casadi::MX &expr, MXPlaceholder::evaluation_mode mode)
            {
                return method->fill_placeholder(type, expr, mode);
            }
            casadi::MX new_placeholder_expression(const casadi::MX &expr, PlaceHolderType type, const std::string &name)
            {
                auto placeholder = casadi::MX::sym(name, expr.size1(), expr.size2());
                placeholders[placeholder] = MXPlaceholder(expr, type, this);
                return placeholder;
            };
            void transcribe(const int K)
            {
                method->transcribe(K);
            }
            StageProblem *parent;
            StageProblem *child;
            Placeholders placeholders;
            std::shared_ptr<StageMethod> method;
            std::map<casadi::MX, casadi::MX, comp_mx> x_next;
            StageProblemDimensions dims;
            std::vector<StageMX> constraints;
            std::vector<StageMX> objective_terms;
            std::vector<casadi::MX> states;
            std::vector<casadi::MX> controls;
        };
    }
} // namespace fatrop_specification
