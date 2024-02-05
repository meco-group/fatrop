#pragma once
#include <casadi/casadi.hpp>
#include <casadi/core/function_internal.hpp>
#include <string>
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
#include <functional>
namespace fatrop
{
    namespace spectool
    {
        typedef int (*eval_t)(const double **arg, double **res,
                              long long int *iw, double *w, int);
        namespace cs = casadi;
        struct CasadiFEJit
        {
            CasadiFEJit(const cs::Function &F, const cs::Dict &jit_options_);
            ~CasadiFEJit();
            casadi::Importer compiler_;
            std::string jit_name_;
            casadi::Dict jit_options_;
            std::string jit_directory;
            bool compiled_jit = false;
            eval_t eval_;
        };

    } // namespace spectool
} // namespace fatrop