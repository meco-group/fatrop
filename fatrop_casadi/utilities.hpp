#pragma once

#include <casadi/casadi.hpp>
namespace fatrop
{
    namespace fatrop_casadi
    {
        struct comp_mx
        {
            bool operator()(const casadi::MX &a, const casadi::MX &b) const
            {
                return a.get() < b.get();
            }
        };
    }
}