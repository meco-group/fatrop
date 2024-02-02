#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
namespace fatrop
{
    namespace spectool
    {
        namespace cs = casadi;
        struct Jacobian
        {
            cs::MX Jx;
            cs::MX x;
            bool is_empty() const {return Jx.is_empty();}
        };
    } // namespace spectrop
} // namespace fatrop