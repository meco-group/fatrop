#pragma once
#include <casadi/casadi.hpp>
#include <vector>
// #include <utility>
#include <functional>
#include <string>
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
        struct StageMX : public casadi::MX
        {
            StageMX(const casadi::MX &expr, bool at_t0, bool at_path, bool at_tf) : casadi::MX(expr), at_t0(at_t0), at_path(at_path), at_tf(at_tf){};
            bool at_t0;
            bool at_path;
            bool at_tf;
        };
        enum PlaceHolderType
        {
            at_t0,
            at_path,
            at_tf,
        };
        class StageProblem;
        struct MXPlaceholder : public casadi::MX
        {
            MXPlaceholder() : casadi::MX(){};
            MXPlaceholder(const casadi::MX &expr, PlaceHolderType type, StageProblem *stage) : casadi::MX(expr), type(type), stage(stage){};
            PlaceHolderType type;
            StageProblem *stage; // TODO check if this can be done in a safer way
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
                throw std::runtime_error("No method set");
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
            casadi::MX operator()(const casadi::MX &expr, MXPlaceholder::evaluation_mode mode);

        public:
            std::vector<std::pair<MXPlaceholder, casadi::MX>> get_all_placeholders(const casadi::MX &expr);
            bool has_placeholders(const casadi::MX &expr, const std::vector<PlaceHolderType> &types);
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

casadi::MX Placeholders::operator()(const casadi::MX &expr, MXPlaceholder::evaluation_mode mode)
{
    casadi::MX ret = expr;
    while (!get_all_placeholders(ret).empty()) // todo re-use this result
    {
        std::vector<casadi::MX> from;
        std::vector<casadi::MX> to;
        auto placeholders = get_all_placeholders(ret);
        while (!placeholders.empty())
        {
            auto p = placeholders.back();
            placeholders.pop_back();
            from.push_back(p.second);
            to.push_back(p.first.stage->fill_placeholder(p.first.type, p.first, mode));
            // std::cout << "replacing " << p.second << " with " << to.back() << std::endl;
        }
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
    // std::cout << "number of placeholders: " << ret.size() << std::endl;
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