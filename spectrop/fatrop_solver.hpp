#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "casadi_utilities.hpp"
#include "solver.hpp"
#include "fatrop_function.hpp"
#include "fatrop_ocp_impl.hpp"
namespace fatrop
{
    namespace spectrop
    {
        namespace cs = casadi;
        struct StageEvaluator
        {
        };
        class SolverFatrop : public SolverInterface
        {
        public:
            void transcribe(const Ocp &ocp_) 
            {
                fatrop_impl = std::make_shared<FatropOcpImpl>(ocp_);
            }
            FatropFunction to_function(const Ocp &ocp_, std::vector<cs::MX> &gist_in, std::vector<cs::MX> &gist_out)
            {
                return FatropFunction(fatrop_impl);
            };
            std::shared_ptr<FatropOcpImpl> fatrop_impl;
        };
    } // namespace spectrop
} // namespace fatrop