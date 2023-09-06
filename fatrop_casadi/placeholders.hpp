#pragma once
#include <casadi/casadi.hpp>
#include <vector>
// #include <utility>
#include <functional>
#include <string>
#include <map>
#include <cassert>
#include "utilities.hpp"
namespace fatrop
{
    namespace fatrop_casadi
    {
        class StageProblem;
        struct StageMX : public casadi::MX
        {
            StageMX(const casadi::MX &expr, bool at_t0, bool at_path, bool at_tf) : casadi::MX(expr), at_t0(at_t0), at_path(at_path), at_tf(at_tf){};
            bool at_t0;
            bool at_path;
            bool at_tf;
        };
        enum class PlaceHolderType
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
        class Placeholders : public std::map<casadi::MX, MXPlaceholder, comp_mx>
        {
        public:
            std::vector<PlaceHolderType> get_all_types(const casadi::MX &expr);
            casadi::MX operator()(const casadi::MX &expr, MXPlaceholder::evaluation_mode mode);

        public:
            std::vector<std::pair<MXPlaceholder, casadi::MX>> get_all_placeholders(const casadi::MX &expr);
            bool has_placeholders(const casadi::MX &expr, const std::vector<PlaceHolderType> &types);
        };
    } // namespace fatrop_casadi
} // namespace fatrop
