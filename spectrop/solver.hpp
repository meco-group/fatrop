
#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "casadi_utilities.hpp"
namespace fatrop
{
    namespace spectrop
    {
        namespace cs = casadi;
        class Ocp;
        class SolverInterface
        {
        public:
            virtual void transcribe(const Ocp &ocp_) = 0;
        };
    } // namespace spectrop
} // namespace fatrop