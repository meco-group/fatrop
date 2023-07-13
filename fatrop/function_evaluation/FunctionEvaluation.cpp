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
#include "function_evaluation/FunctionEvaluation.hpp"
using namespace fatrop;
using namespace std;

fatrop_int EvalBase::eval_bf(const double **arg, MAT *bf_mat)
{
#if DEBUG
    assert(bf_mat->m >= out_m);
    assert(bf_mat->n >= out_n);
#endif
    #ifndef ENABLE_MULTITHREADING
    double *buffer_p = buffer.data();
    #else
    double *buffer_p = buffer[omp_get_thread_num()].data();
    #endif
    // todo make this static polymorphism using CRTP
    fatrop_int res = eval_buffer(arg);
    PACKMAT(out_m, out_n, buffer_p, out_m, bf_mat, 0, 0);
    return res;
}
fatrop_int EvalBase::eval_array(const double **arg, double *array)
{
    #ifndef ENABLE_MULTITHREADING
    double *buffer_p = buffer.data();
    #else
    double *buffer_p = buffer[omp_get_thread_num()].data();
    #endif
    // todo make this static polymorphism using CRTP
    fatrop_int res = eval_buffer(arg);
    memcpy(array, buffer_p, out_nnz * sizeof(double));
    return res;
}