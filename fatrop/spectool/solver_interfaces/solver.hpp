
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
            virtual cs::Dict stats() = 0;
            virtual void transcribe(const Ocp &ocp_, const cs::Dict &opts) = 0;
            virtual cs::Function to_function(const std::string &name, const Ocp &ocp, std::vector<cs::MX> &gist_in, std::vector<cs::MX> &gist_out, const cs::Dict &opts) = 0;
            virtual void gist(const Ocp &ocp_, std::vector<cs::MX> &in, std::vector<cs::MX> &out) = 0;
            // make the destructor virtual, so that the derived class's destructor is called
            virtual ~SolverInterface() {}
        };
    } // namespace spectrop
} // namespace fatrop