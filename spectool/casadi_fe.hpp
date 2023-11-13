#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "casadi_utilities.hpp"
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "casadi_jit.hpp"
namespace fatrop
{
    namespace spectool
    {
        namespace cs = casadi;
        class CasadiFEWrap
        {
        public:
            CasadiFEWrap(){};
            CasadiFEWrap(const cs::Function &func, bool expand, bool jit, const cs::Dict& jit_options_) : func_(expand ? func.expand() : func), jit_(jit)
            {
                m = (int)func_.size1_out(0);
                n = (int)func_.size2_out(0);
                // func = func.expand();
                mem = 0;
                // allocate work vectors
                size_t sz_arg,
                    sz_res, sz_iw, sz_w;
                sz_arg = func.n_in();
                sz_res = func.n_out();
                func.sz_work(sz_arg, sz_res, sz_iw, sz_w);
                iw.resize(sz_iw);
                w.resize(sz_w);
                bufout.resize(func.nnz_out(0));
                bufdata.resize(sz_res > 0 ? sz_res : 1);
                resdata.resize(sz_res > 0 ? sz_res : 1);
                argdata.resize(sz_arg > 0 ? sz_arg : 1);
                n_in = func.n_in();
                if (jit)
                {
                    fe_jit_ = std::make_unique<CasadiFEJit>(func_, jit_options_);
                }
                // // assert dense matrix output
                // assert(func.nnz_out(0) == m * n);
            };
            void eval(const double **args, MAT *res)
            {
                // assert(res->m == m);
                // assert(res->n == n);
                eval(args, bufout.data());
                PACKMAT(m, n, bufout.data(), m, res, 0, 0);
            };
            void eval(const double **args, double *res)
            {
                // inputs
                for (int j = 0; j < n_in; j++)
                    argdata[j] = args[j];
                // outputs
                if (res)
                {
                    resdata[0] = res;
                    jit_ ? fe_jit_->eval_(argdata.data(), resdata.data(), iw.data(), w.data(), 0) : func_(argdata.data(), resdata.data(), iw.data(), w.data(), 0);
                }
            };

        private:
            cs::Function func_;
            int m;
            int n;
            int mem;
            int n_in;
            bool jit_;
            std::vector<double> bufout;
            std::vector<double *> bufdata;
            std::vector<double *> resdata;
            std::vector<const double *> argdata;
            std::vector<long long int> iw;
            std::vector<double> w;
            std::unique_ptr<CasadiFEJit> fe_jit_;
        };
        // implementation of OCPAbstract, given an OCP
    } // namespace spectool
} // namespace fatrop