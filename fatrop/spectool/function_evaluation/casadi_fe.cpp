#include "casadi_fe.hpp"

namespace fatrop
{
    namespace spectool
    {
        CasadiFEWrap::CasadiFEWrap(const cs::Function &func, bool expand, bool jit, const cs::Dict &jit_options_, CasadiJitCache &cache_map) : func_(expand ? func.expand() : func), jit_(jit)
        {
            m = (int)func_.size1_out(0);
            n = (int)func_.size2_out(0);
            // func = func.expand();
            mem = 0;
            // allocate work vectors
            size_t sz_arg,
                sz_res, sz_iw, sz_w;
            sz_arg = func_.n_in();
            sz_res = func_.n_out();
            func_.sz_work(sz_arg, sz_res, sz_iw, sz_w);
#ifndef ENABLE_MULTITHREADING
            iw.resize(sz_iw);
            w.resize(sz_w);
            bufout.resize(func_.nnz_out(0));
            bufdata.resize(sz_res > 0 ? sz_res : 1);
            resdata.resize(sz_res > 0 ? sz_res : 1);
            argdata.resize(sz_arg > 0 ? sz_arg : 1);
#else
            auto num_threads = omp_get_max_threads();
            for (int thread_i = 0; thread_i < num_threads; thread_i++)
            {
                iw[thread_i].resize(sz_iw);
                w[thread_i].resize(sz_w);
                bufout[thread_i].resize(func_.nnz_out(0));
                bufdata[thread_i].resize(sz_res > 0 ? sz_res : 1);
                resdata[thread_i].resize(sz_res > 0 ? sz_res : 1);
                argdata[thread_i].resize(sz_arg > 0 ? sz_arg : 1);
            }
#endif

            n_in = func_.n_in();
            if (jit)
            {
                // compute the function hash
                size_t func_hash = std::hash<std::string>{}(func_.serialize());
                // check if the function is available in the cache
                if (cache_map.find(func_hash) != cache_map.end())
                {
                    fe_jit_ = cache_map[func_hash];
                }
                else
                {
                    fe_jit_ = std::make_shared<CasadiFEJit>(func_, jit_options_);
                    // add the function to the cache
                    cache_map[func_hash] = fe_jit_;
                }
            }
            // // assert dense matrix output
            // assert(func.nnz_out(0) == m * n);
        };
        void CasadiFEWrap::eval(const double **args, MAT *res)
        {
// assert(res->m == m);
// assert(res->n == n);
#ifndef ENABLE_MULTITHREADING
            eval(args, bufout.data());
            PACKMAT(m, n, bufout.data(), m, res, 0, 0);
#else
            int tid = omp_get_thread_num();
            eval(args, bufout[tid].data());
            PACKMAT(m, n, bufout[tid].data(), m, res, 0, 0);
#endif
        };
        void CasadiFEWrap::eval(const double **args, double *res)
        {
#ifndef ENABLE_MULTITHREADING
            // inputs
            for (int j = 0; j < n_in; j++)
                argdata[j] = args[j];
            // outputs
            if (res)
            {
                resdata[0] = res;
                jit_ ? fe_jit_->eval_(argdata.data(), resdata.data(), iw.data(), w.data(), 0) : func_(argdata.data(), resdata.data(), iw.data(), w.data(), 0);
            }
#else
            int tid = omp_get_thread_num();
            // inputs
            for (int j = 0; j < n_in; j++)
                argdata[tid][j] = args[j];
            // outputs
            if (res)
            {
                resdata[tid][0] = res;
                jit_ ? fe_jit_->eval_(argdata[tid].data(), resdata[tid].data(), iw[tid].data(), w[tid].data(), 0) : func_(argdata[tid].data(), resdata[tid].data(), iw[tid].data(), w[tid].data(), 0);
            }
#endif
        };

    }
}