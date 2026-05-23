/*
 *    MIT No Attribution
 *
 *    Copyright 2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
 *
 *    Permission is hereby granted, free of charge, to any person obtaining a copy of this
 *    software and associated documentation files (the "Software"), to deal in the Software
 *    without restriction, including without limitation the rights to use, copy, modify,
 *    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 *    permit persons to whom the Software is furnished to do so.
 *
 *    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 *    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 *    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 *    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

/**
 *  Example program demonstrating the usage of C-code generated from CasADi
 *  Generated code is encapsulated in self-contained code with the entry points
 *  defined below.
 *  Note that other usage, e.g. accessing the internal data structures in the
 *  generated files is not recommended and subject to change.
 *
 *  We show how the generated code can be used from C (or C++), without requiring
 *  any CasADi classes as well as how to use it from C++ using CasADi's external.
 *
 *  Joel Andersson, 2013-2015
 */

#include <cassert>
#include <stdio.h>
#include <string>
#include <vector>
#include "casadi_generated.h"

int main()
{
    printf("---\n");
    printf("Standalone usage from C/C++:\n");
    printf("\n");

    // Get number of inputs and outputs
    casadi_int n_in = opti_func_n_in();
    casadi_int n_out = opti_func_n_out();

    // Get sizes of the required work vectors
    casadi_int sz_arg = n_in, sz_res = n_out, sz_iw = 0, sz_w = 0;
    opti_func_work(&sz_arg, &sz_res, &sz_iw, &sz_w);
    printf("Work vector sizes:\n");
    printf("sz_arg = %lld, sz_res = %lld, sz_iw = %lld, sz_w = %lld\n\n", sz_arg, sz_res, sz_iw, sz_w);

    // Allocate input/output buffers and work vectors
    const double *arg[sz_arg];
    double *res[sz_res];
    casadi_int iw[sz_iw];
    double w[sz_w];

    // Function input and output
    std::vector<double> res_buffer(opti_func_sparsity_out(0)[0]);

    // Set output buffer
    res[0] = res_buffer.data();

    // Allocate memory (thread-safe)
    opti_func_incref();

    // Checkout thread-local memory (not thread-safe)
    int mem = opti_func_checkout();

    // Evaluation is thread-safe
    if (opti_func(arg, res, iw, w, mem))
    {
        printf("Error during function evaluation\n");
        return 1;
    }

    // Release thread-local (not thread-safe)
    opti_func_release(mem);

    // Free memory (thread-safe)
    opti_func_decref();

    printf("Function evaluation completed successfully\n");
    return 0;
}
