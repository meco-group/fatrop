#pragma once
#include <casadi/casadi.hpp>
#include <vector>
// #include <utility>
#include <functional>
#include <map>
#include "fatrop-casadi-problem.hpp"
#include "utilities.hpp"
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
            path,
            path_t0tf,
            path_t0,
            path_tf,
            sum,      // include_first = False, inlclude_last = False
            sum_t0tf, // include_first = True, include_last = True
            sum_t0,   // include_first = True, include_last = False
            sum_tf,   // include_first = False, include_last = True
        };
        class StageProblemInternal;
        struct MXPlaceholder : public casadi::MX
        {
            MXPlaceholder() : casadi::MX(){};
            MXPlaceholder(const casadi::MX &expr, PlaceHolderType type, StageProblemInternal *stage) : casadi::MX(expr), type(type), stage(stage){};
            PlaceHolderType type;
            StageProblemInternal *stage; // TODO check if this can be done in a safer way
            enum evaluation_mode
            {
                transcribe,
                evaluate
            };
        };
        class StageMethod
        {
        public:
            virtual casadi::MX fill_placeholder(const PlaceHolderType type, const casadi::MX &expr, MXPlaceholder::evaluation_mode mode)
            {
                std::runtime_error("No method set");
                return casadi::MX();
            };
            virtual void transcribe(const int K)
            {
                std::runtime_error("No method set");
            };
        };

        class Placeholders : public std::map<casadi::MX, MXPlaceholder, comp_mx>
        {
        public:
            std::vector<PlaceHolderType> get_all_types(const casadi::MX &expr);
            casadi::MX operator()(casadi::MX &expr, MXPlaceholder::evaluation_mode mode);

        public:
            std::vector<std::pair<MXPlaceholder, casadi::MX>> get_all_placeholders(const casadi::MX &expr);
            bool has_placeholders(const casadi::MX &expr, const std::vector<PlaceHolderType> &types);
        };

        struct StageProblem : public SharedObj<StageProblemInternal, StageProblem>
        {
            using SharedObj<StageProblemInternal, StageProblem>::SharedObj;
        };

        class StageProblemInternal
        {
        public:
            casadi::MX state(const int m, const int n = 1)
            {
                casadi::MX ret = casadi::MX::sym("x", m, n);
                states.push_back(ret);
                // x_next[ret] = casadi::MX::sym("dummy", 0, 0);
                return ret;
            }
            casadi::MX control(const int m, const int n = 1)
            {
                casadi::MX ret = casadi::MX::sym("u", m, n);
                controls.push_back(ret);
                return ret;
            }
            void set_next(const casadi::MX &x, const casadi::MX &x_next)
            {
                this->x_next[x] = x_next;
            }
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
                    return new_placeholder_expression(expr, PlaceHolderType::sum_t0tf);
                else if (include_first && !include_last)
                    return new_placeholder_expression(expr, PlaceHolderType::sum_t0);
                else if (!include_first && include_last)
                    return new_placeholder_expression(expr, PlaceHolderType::sum_tf);
                else if (!include_first && !include_last)
                    return new_placeholder_expression(expr, PlaceHolderType::sum);
            }
            void subject_to(const casadi::MX &expr, bool include_first = false, bool include_last = false)
            {
                if (placeholders.has_placeholders(expr, {PlaceHolderType::at_t0, PlaceHolderType::at_tf})) // check for point constraints
                {
                    constraints.push_back(new_placeholder_expression(expr, PlaceHolderType::path));
                    return;
                }
                if (include_first && include_last)
                    constraints.push_back(new_placeholder_expression(expr, PlaceHolderType::path_t0tf));
                else if (include_first && !include_last)
                    constraints.push_back(new_placeholder_expression(expr, PlaceHolderType::path_t0));
                else if (!include_first && include_last)
                    constraints.push_back(new_placeholder_expression(expr, PlaceHolderType::path_tf));
                else if (!include_first && !include_last)
                    constraints.push_back(new_placeholder_expression(expr, PlaceHolderType::path));
            }
            void add_objective(const casadi::MX &expr)
            {
                objective_terms.push_back(expr);
            }
            casadi::MX fill_placeholder(const PlaceHolderType type, const casadi::MX &expr, MXPlaceholder::evaluation_mode mode)
            {
                return method->fill_placeholder(type, expr, mode);
            }
            casadi::MX new_placeholder_expression(const casadi::MX &expr, PlaceHolderType type)
            {
                casadi::MX placeholder;
                placeholders[placeholder] = MXPlaceholder(expr, type, this);
                return placeholder;
            };
            void transcribe(const int K)
            {
                method->transcribe(K);
            }
            StageProblem parent;
            StageProblem child;
            Placeholders placeholders;
            std::shared_ptr<StageMethod> method;
            std::map<casadi::MX, casadi::MX, comp_mx> x_next;
            StageProblemDimensions dims;
            std::vector<casadi::MX> constraints;
            std::vector<casadi::MX> objective_terms;
            std::vector<casadi::MX> states;
            std::vector<casadi::MX> controls;
        };
    }
} // namespace fatrop_specification

///////////// IMPLEMENTATION
using namespace fatrop::fatrop_casadi;

std::vector<PlaceHolderType> Placeholders::get_all_types(const casadi::MX &expr)
{
    std::vector<PlaceHolderType> ret;
    for (const auto &p : get_all_placeholders(expr))
    {
        ret.push_back(p.first.type);
    }
    return ret;
}

casadi::MX Placeholders::operator()(casadi::MX &expr, MXPlaceholder::evaluation_mode mode)
{
    casadi::MX ret = expr;
    while (!get_all_placeholders(ret).empty()) // todo re-use this result
    {
        std::vector<casadi::MX> from;
        std::vector<casadi::MX> to;
        auto placeholders = get_all_placeholders(ret);
        auto p = placeholders.back();
        placeholders.pop_back();
        from.push_back(p.second);
        to.push_back(p.first.stage->fill_placeholder(p.first.type, p.second, mode));
        ret = casadi::MX::substitute({ret}, from, to)[0];
    }
    return ret;
}
std::vector<std::pair<MXPlaceholder, casadi::MX>> Placeholders::get_all_placeholders(const casadi::MX &expr)
{
    auto ret = std::vector<std::pair<MXPlaceholder, casadi::MX>>();
    auto syms = casadi::MX::symvar(expr);
    for (auto &sym : syms)
    {
        if (find(sym) != end())
        {
            ret.push_back(std::make_pair(operator[](sym), sym));
        }
    }
    return ret;
}
bool Placeholders::has_placeholders(const casadi::MX &expr, const std::vector<PlaceHolderType> &types)
{
    auto placeholders = get_all_placeholders(expr);
    for (auto &p : placeholders)
    {
        if (std::find(types.begin(), types.end(), p.first.type) != types.end())
            return true;
    }
    return false;
}