#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
namespace fatrop
{
    namespace spectool
    {
        namespace cs = casadi;
        struct Hessian
        {
            cs::MX Hxx;
            cs::MX Hx;
            cs::MX x;
            cs::MX lam;
            bool is_empty() const {return Hxx.is_empty();}
        };
    } // namespace spectrop
} // namespace fatrop