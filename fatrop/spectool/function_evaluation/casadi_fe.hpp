#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
#include "fatrop/blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "casadi_jit.hpp"
#include "casadi_jit_cache.hpp"
#ifdef ENABLE_MULTITHREADING
#include <omp.h>
#endif
namespace fatrop
{
    namespace spectool
    {
        namespace cs = casadi;
        class CasadiFEWrap
        {
        public:
            CasadiFEWrap(){};
            CasadiFEWrap(const cs::Function &func, bool expand, bool jit, const cs::Dict &jit_options_, CasadiJitCache &cache_map);
            void eval(const double **args, MAT *res);
            void eval(const double **args, double *res);

        private:
            cs::Function func_;
            int m;
            int n;
            int mem;
            int n_in;
            bool jit_;
#ifndef ENABLE_MULTITHREADING
            std::vector<double> bufout;
            std::vector<double *> bufdata;
            std::vector<double *> resdata;
            std::vector<const double *> argdata;
            std::vector<long long int> iw;
            std::vector<double> w;
#else
            std::vector<std::vector<double>>  bufout = std::vector<std::vector<double>>(omp_get_max_threads()) ;
            std::vector<std::vector<double *>>  bufdata = std::vector<std::vector<double *>>(omp_get_max_threads()) ;
            std::vector<std::vector<double *>>  resdata = std::vector<std::vector<double *>>(omp_get_max_threads()) ;
            std::vector<std::vector<const double *>>  argdata = std::vector<std::vector<const double *>>(omp_get_max_threads()) ;
            std::vector<std::vector<long long int>>  iw = std::vector<std::vector<long long int>>(omp_get_max_threads()) ;
            std::vector<std::vector<double>>  w = std::vector<std::vector<double>>(omp_get_max_threads()) ;
#endif
            std::shared_ptr<CasadiFEJit> fe_jit_;
        };
        // implementation of OCPAbstract, given an OCP
    } // namespace spectool
} // namespace fatrop