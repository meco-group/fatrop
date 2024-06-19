/**
*  THIS FILE IS MODIFILED FROM blasfeo_ref/x_blas_ref3.c from the BLASFEO library.

WITH THE FOLLOWING COPYRIGHT NOTICE:



BLASFEO -- BLAS For Embedded Optimization.
Copyright (C) 2019 by Gianluca Frison.
Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.
All rights reserved.

The 2-Clause BSD License

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**/

extern "C"
{
#include <blasfeo.h>
}
#include <iostream>

#if defined(BLASFEO_REF_API)
#else
#define MATEL BLASFEO_DMATEL
#define MAT blasfeo_dmat
#define VEC blasfeo_dvec
void blasfeo_ref_dtrsm_rlnn_copy(int m, int n, double alpha, struct MAT *sA, int ai, int aj, struct MAT *sB, int bi, int bj, struct MAT *sD, int di, int dj)
{
    if (m <= 0 | n <= 0)
        return;

    // invalidate stored inverse diagonal of result matrix
    sD->use_dA = 0;

    char cl = 'l';
    char cn = 'n';
    char cr = 'r';
    char ct = 't';
    char cu = 'u';
    int i1 = 1;
    int ii, jj, kk, id;
    double
        d_00,
        d_01,
        d_10, d_11;
    int aai = ai;
    int aaj = aj;
    int bbi = bi;
    int bbj = bj;
    int ddi = di;
    int ddj = dj;
    double *dA = sA->dA;
    if (ai == 0 & aj == 0)
    {
        if (sA->use_dA < n)
        {
            // invert diagonal of pA
            for (ii = 0; ii < n; ii++)
                dA[ii] = 1.0 / MATEL(sA, aai + ii, aaj + ii);
            // use only now
            sA->use_dA = n;
        }
    }
    else
    {
        for (ii = 0; ii < n; ii++)
            dA[ii] = 1.0 / MATEL(sA, aai + ii, aaj + ii);
        sA->use_dA = 0; // nonzero offset makes diagonal dirty
    }
    // solve
    jj = 0;
    for (; jj < n - 1; jj += 2)
    {
        ii = 0;
        id = n - jj - 2;
        for (; ii < m - 1; ii += 2)
        {
            d_00 = alpha * MATEL(sB, bbi + ii + 0, bbj + (id + 0));
            d_10 = alpha * MATEL(sB, bbi + ii + 1, bbj + (id + 0));
            d_01 = alpha * MATEL(sB, bbi + ii + 0, bbj + (id + 1));
            d_11 = alpha * MATEL(sB, bbi + ii + 1, bbj + (id + 1));
            kk = id + 2;
            for (; kk < n; kk++)
            {
                d_00 -= MATEL(sA, aai + kk + 0, aaj + (id + 0)) * MATEL(sD, ddi + ii + 0, ddj + (kk + 0));
                d_10 -= MATEL(sA, aai + kk + 0, aaj + (id + 0)) * MATEL(sD, ddi + ii + 1, ddj + (kk + 0));
                d_01 -= MATEL(sA, aai + kk + 0, aaj + (id + 1)) * MATEL(sD, ddi + ii + 0, ddj + (kk + 0));
                d_11 -= MATEL(sA, aai + kk + 0, aaj + (id + 1)) * MATEL(sD, ddi + ii + 1, ddj + (kk + 0));
            }
            d_01 *= dA[id + 1];
            d_11 *= dA[id + 1];
            d_00 -= MATEL(sA, aai + id + 1, aaj + (id + 0)) * d_01;
            d_10 -= MATEL(sA, aai + id + 1, aaj + (id + 0)) * d_11;
            d_00 *= dA[id + 0];
            d_10 *= dA[id + 0];
            MATEL(sD, ddi + ii + 0, ddj + (id + 0)) = d_00;
            MATEL(sD, ddi + ii + 1, ddj + (id + 0)) = d_10;
            MATEL(sD, ddi + ii + 0, ddj + (id + 1)) = d_01;
            MATEL(sD, ddi + ii + 1, ddj + (id + 1)) = d_11;
        }
        for (; ii < m; ii++)
        {
            d_00 = alpha * MATEL(sB, bbi + ii + 0, bbj + (id + 0));
            d_01 = alpha * MATEL(sB, bbi + ii + 0, bbj + (id + 1));
            kk = id + 2;
            for (; kk < n; kk++)
            {
                d_00 -= MATEL(sA, aai + kk + 0, aaj + (id + 0)) * MATEL(sD, ddi + ii + 0, ddj + (kk + 0));
                d_01 -= MATEL(sA, aai + kk + 0, aaj + (id + 1)) * MATEL(sD, ddi + ii + 0, ddj + (kk + 0));
            }
            d_01 *= dA[id + 1];
            d_00 -= MATEL(sA, aai + id + 1, aaj + (id + 0)) * d_01;
            d_00 *= dA[id + 0];
            MATEL(sD, ddi + ii + 0, ddj + (id + 0)) = d_00;
            MATEL(sD, ddi + ii + 0, ddj + (id + 1)) = d_01;
        }
    }
    for (; jj < n; jj++)
    {
        ii = 0;
        id = n - jj - 1;
        for (; ii < m - 1; ii += 2)
        {
            d_00 = alpha * MATEL(sB, bbi + ii + 0, bbj + (id + 0));
            d_10 = alpha * MATEL(sB, bbi + ii + 1, bbj + (id + 0));
            kk = id + 1;
            for (; kk < n; kk++)
            {
                d_00 -= MATEL(sA, aai + kk + 0, aaj + (id + 0)) * MATEL(sD, ddi + ii + 0, ddj + (kk + 0));
                d_10 -= MATEL(sA, aai + kk + 0, aaj + (id + 0)) * MATEL(sD, ddi + ii + 1, ddj + (kk + 0));
            }
            d_00 *= dA[id + 0];
            d_10 *= dA[id + 0];
            MATEL(sD, ddi + ii + 0, ddj + (id + 0)) = d_00;
            MATEL(sD, ddi + ii + 1, ddj + (id + 0)) = d_10;
        }
        for (; ii < m; ii++)
        {
            d_00 = alpha * MATEL(sB, bbi + ii, bbj + (id));
            kk = id + 1;
            for (; kk < n; kk++)
                d_00 -= MATEL(sA, aai + kk, aaj + (id)) * MATEL(sD, ddi + ii, ddj + (kk));
            MATEL(sD, ddi + ii, ddj + (id)) = d_00 * dA[id];
        }
    }
    return;
}
#endif