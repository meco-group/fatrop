#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <map>
namespace fatrop
{
    namespace spectool
    {
        struct MxHash
        {
            std::size_t operator()(const casadi::MX &mx) const
            {
                return std::hash<casadi::MXNode*>()(mx.get());
            }
        };
        struct MxEqual
        {
            bool operator()(const casadi::MX &mx1, const casadi::MX &mx2) const
            {
                return mx1.get() == mx2.get();
            }
        };
        typedef std::unordered_set<casadi::MX, MxHash, MxEqual> uo_set_mx;
        // typedef std::set<casadi::MX, MxHash, MxEqual> set_mx;
        template <typename T>
        struct uo_map_mx: public std::unordered_map<casadi::MX, T, MxHash, MxEqual>
        {
        };

        typedef std::unordered_map<casadi::MX, casadi::MX, MxHash, MxEqual> uo_map_mx_mx;
        // template <typename T>
        // struct map_mx: public std::map<casadi::MX, T, MxHash, MxEqual>
        // {
        // };
    } // namespace spectool
} // namespace fatrop