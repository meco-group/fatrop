
#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
namespace fatrop
{
    namespace spectool
    {
        namespace cs = casadi;
        class Ocp;
        class SolverInterface
        {
        public:
            // virtual void transcribe(const Ocp &ocp_) = 0;
            // virtual cs::Function to_function(const Ocp &ocp_, std::vector<cs::MX>& gist_in,std::vector<cs::MX>& gist_out) = 0;
        };
    } // namespace spectrop
} // namespace fatrop