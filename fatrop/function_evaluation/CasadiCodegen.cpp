/*
 * Fatrop - A fast trajectory optimization solver
 * Copyright (C) 2022, 2023 Lander Vanroye, KU Leuven. All rights reserved.
 *
 * This file is part of Fatrop.
 *
 * Fatrop is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fatrop is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Fatrop.  If not, see <http://www.gnu.org/licenses/>. */
#include "function_evaluation/CasadiCodegen.hpp"
using namespace fatrop;
using namespace std;

EvalCasGen::EvalCasGen(const shared_ptr<DLHandler> &handle_, const std::string &function_name) : handle(handle_)
{
    void *handle_p = handle_->handle;

    /* Memory management -- increase reference counter */
    incref = (signal_t)dlsym(handle_p, (function_name + (std::string) "_incref").c_str());
    if (dlerror())
        dlerror(); // No such function, reset error flags

    /* Memory management -- decrease reference counter */
    decref = (signal_t)dlsym(handle_p, (function_name + (std::string) "_decref").c_str());
    if (dlerror())
        dlerror(); // No such function, reset error flags

    /* Thread-local memory management -- checkout memory */
    casadi_checkout_t checkout = (casadi_checkout_t)dlsym(handle_p, (function_name + (std::string) "_checkout").c_str());
    if (dlerror())
        dlerror(); // No such function, reset error flags

    /* T memory management -- release memory */
    release = (casadi_release_t)dlsym(handle_p, (function_name + (std::string) "_release").c_str());
    if (dlerror())
        dlerror(); // No such function, reset error flags

    /* Number of inputs */
    getint_t n_in_fcn = (getint_t)dlsym(handle_p, (function_name + (std::string) "_n_in").c_str());
    // if (dlerror()) return 1;
    n_in = n_in_fcn();

    /* Number of outputs */
    getint_t n_out_fcn = (getint_t)dlsym(handle_p, (function_name + (std::string) "_n_out").c_str());
    // if (dlerror()) return 1;
    casadi_int n_out = n_out_fcn();
    assert(n_out == 1);

    // Checkout thread-local memory (not thread-safe)
    // Note MAX_NUM_THREADS
    mem = checkout();

    /* Get sizes of the required work vectors */
    casadi_int sz_arg = n_in, sz_res = n_out, sz_iw = 0, sz_w = 0;
    work_t work = (work_t)dlsym(handle_p, (function_name + (std::string) "_work").c_str());

    if (dlerror())
        dlerror(); // No such function, reset error flags
    assert((work && work(&sz_arg, &sz_res, &sz_iw, &sz_w)) == 0);
    work_vector_d.resize(sz_w);
    work_vector_i.resize(sz_iw);
    w = work_vector_d.data();
    iw = work_vector_i.data();
    /* Input sparsities */
    // sparsity_t sp_in = (sparsity_t)dlsym(handle, (function_name + (std::string) "_sparsity_in").c_str());
    // assert(dlerror() == 0);

    /* Output sparsities */
    sparsity_t sp_out = (sparsity_t)dlsym(handle_p, (function_name + (std::string) "_sparsity_out").c_str());
    assert(dlerror() == 0);
    const casadi_int *sparsity_out_ci = sp_out(0); // ci stands for casadi_int
    out_m = sparsity_out_ci[0];
    out_n = sparsity_out_ci[1];
    out_nnz = sparsity_out_ci[out_n + 2];
    sparsity_out.resize(2 + out_n + 1 + out_nnz, 0);
    sparsity_out.assign(sparsity_out_ci, sparsity_out_ci + 2 + out_n + 1 + out_nnz);

    /* Function for numerical evaluation */
    eval = (eval_t)dlsym(handle_p, function_name.c_str());
    if (dlerror())
    {
        printf("Failed to retrieve \"f\" function.\n");
    }

// allocate output buffer
#ifndef ENABLE_MULTITHREADING
    buffer.resize(out_nnz, 0.0);
    output_buffer_p = buffer.data();
#else
    for (int i = 0; i < omp_get_max_threads(); i++)
    {
        buffer[i].resize(out_nnz, 0.0);
        output_buffer_p[i] = buffer[i].data();
    }
#endif
}
EvalCasGen::EvalCasGen(
    const signal_t incref,
    const signal_t decref,
    const casadi_checkout_t checkout,
    const casadi_release_t release,
    const getint_t n_in_fcn,
    const getint_t n_out_fcn,
    const sparsity_t sp_in,
    const sparsity_t sp_out,
    const work_t work,
    const eval_t eval):eval(eval), incref(incref), decref(decref), release(release)
{
    /* Memory management -- increase reference counter */

    /* Number of inputs */
    n_in = n_in_fcn();

    /* Number of outputs */
    casadi_int n_out = n_out_fcn();
    assert(n_out == 1);

    // Checkout thread-local memory (not thread-safe)
    // Note MAX_NUM_THREADS
    mem = checkout();

    /* Get sizes of the required work vectors */
    casadi_int sz_arg = n_in, sz_res = n_out, sz_iw = 0, sz_w = 0;

    // assert((work && work(&sz_arg, &sz_res, &sz_iw, &sz_w)) == 0);
    assert(work(&sz_arg, &sz_res, &sz_iw, &sz_w) == 0);
    work_vector_d.resize(sz_w);
    work_vector_i.resize(sz_iw);
    w = work_vector_d.data();
    iw = work_vector_i.data();

    /* Output sparsities */
    const casadi_int *sparsity_out_ci = sp_out(0); // ci stands for casadi_int
    out_m = sparsity_out_ci[0];
    out_n = sparsity_out_ci[1];
    out_nnz = sparsity_out_ci[out_n + 2];
    sparsity_out.resize(2 + out_n + 1 + out_nnz, 0);
    sparsity_out.assign(sparsity_out_ci, sparsity_out_ci + 2 + out_n + 1 + out_nnz);

// allocate output buffer
#ifndef ENABLE_MULTITHREADING
    buffer.resize(out_nnz, 0.0);
    output_buffer_p = buffer.data();
#else
    for (int i = 0; i < omp_get_max_threads(); i++)
    {
        buffer[i].resize(out_nnz, 0.0);
        output_buffer_p[i] = buffer[i].data();
    }
#endif
}

int EvalCasGen::eval_buffer(const double **arg)
{
    w = work_vector_d.data();
    iw = work_vector_i.data();
#ifndef ENABLE_MULTITHREADING
    output_buffer_p = buffer.data();
    if (eval(arg, &output_buffer_p, iw, w, mem))
        return 1;
#else
    int tid = omp_get_thread_num();
    output_buffer_p[tid] = buffer[tid].data();
    // printf("using thread %d\n", tid);
    if (eval(arg, &output_buffer_p[tid], iw, w, mem))
        return 1;
#endif

    return 0;
}

EvalCasGen::~EvalCasGen()
{
    // Release thread-local (not thread-safe)
    release(mem);
    // Free memory (thread-safe)
    decref();
}