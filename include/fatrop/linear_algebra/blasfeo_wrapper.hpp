#ifndef __fatrop_linear_algebra_blasfeo_wrapper_hpp__
#define __fatrop_linear_algebra_blasfeo_wrapper_hpp__

#include "fatrop/common/exception.hpp"
extern "C"
{
#include <blasfeo.h>
}

namespace fatrop
{
    static inline void blasfeo_daxpby_debug(int m, double alpha, blasfeo_dvec *x, int xi,
                                            double beta, blasfeo_dvec *y, int yi, blasfeo_dvec *z,
                                            int zi)
    {
        fatrop_dbg_assert(m > 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && zi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m && yi + m <= y->m && zi + m <= z->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        blasfeo_daxpby(m, alpha, x, xi, beta, y, yi, z, zi);
    }

    static inline void blasfeo_daxpy_debug(int m, double alpha, blasfeo_dvec *x, int xi,
                                           blasfeo_dvec *y, int yi, blasfeo_dvec *z, int zi)
    {
        fatrop_dbg_assert(m > 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && zi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m && yi + m <= y->m && zi + m <= z->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        blasfeo_daxpy(m, alpha, x, xi, y, yi, z, zi);
    }

    static inline void blasfeo_dvecse_debug(int m, double alpha, blasfeo_dvec *x, int xi)
    {
        fatrop_dbg_assert(m > 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && "Vector index must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m &&
                          "Vector index plus size must not exceed vector dimension");
        blasfeo_dvecse(m, alpha, x, xi);
    }

    static inline void blasfeo_dvecsc_debug(int m, double alpha, blasfeo_dvec *x, int xi)
    {
        fatrop_dbg_assert(m > 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && "Vector index must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m &&
                          "Vector index plus size must not exceed vector dimension");
        blasfeo_dvecsc(m, alpha, x, xi);
    }

    static inline void blasfeo_dveccpsc_debug(int m, double alpha, blasfeo_dvec *x, int xi,
                                              blasfeo_dvec *y, int yi)
    {
        fatrop_dbg_assert(m > 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m && yi + m <= y->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        blasfeo_dveccpsc(m, alpha, x, xi, y, yi);
    }

    static inline void blasfeo_dveccp_debug(int m, blasfeo_dvec *x, int xi, blasfeo_dvec *y, int yi)
    {
        fatrop_dbg_assert(m > 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m && yi + m <= y->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        blasfeo_dveccp(m, x, xi, y, yi);
    }

    static inline void blasfeo_dvecmul_debug(int m, blasfeo_dvec *x, int xi, blasfeo_dvec *y,
                                             int yi, blasfeo_dvec *z, int zi)
    {
        fatrop_dbg_assert(m > 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && zi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m && yi + m <= y->m && zi + m <= z->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        blasfeo_dvecmul(m, x, xi, y, yi, z, zi);
    }

    static inline void blasfeo_drowpe_debug(int kmax, int *ipiv,  blasfeo_dmat *sA)
    {
        fatrop_dbg_assert(kmax > 0 && "kmax must be positive");
        blasfeo_drowpe(kmax, ipiv, sA);
    }

    static inline void blasfeo_dvecpe_debug(int kmax, int *ipiv,  blasfeo_dvec *sx, int xi)
    {
        fatrop_dbg_assert(kmax > 0 && "kmax must be positive");
        fatrop_dbg_assert(xi >= 0 && "xi must be non-negative");
        fatrop_dbg_assert(xi + kmax <= sx->m && "xi + kmax must not exceed vector dimension");
        blasfeo_dvecpe(kmax, ipiv, sx, xi);
    }

    static inline void blasfeo_dvecpei_debug(int kmax, int *ipiv,  blasfeo_dvec *sx, int xi)
    {
        fatrop_dbg_assert(kmax > 0 && "kmax must be positive");
        fatrop_dbg_assert(xi >= 0 && "xi must be non-negative");
        fatrop_dbg_assert(xi + kmax <= sx->m && "xi + kmax must not exceed vector dimension");
        blasfeo_dvecpei(kmax, ipiv, sx, xi);
    }

    static inline void blasfeo_drowpei_debug(int kmax, int *ipiv,  blasfeo_dmat *sA)
    {
        fatrop_dbg_assert(kmax > 0 && "kmax must be positive");
        blasfeo_drowpei(kmax, ipiv, sA);
    }

    static inline void blasfeo_dcolpe_debug(int kmax, int *ipiv,  blasfeo_dmat *sA)
    {
        fatrop_dbg_assert(kmax > 0 && "kmax must be positive");
        blasfeo_dcolpe(kmax, ipiv, sA);
    }

    static inline void blasfeo_dcolpei_debug(int kmax, int *ipiv,  blasfeo_dmat *sA)
    {
        fatrop_dbg_assert(kmax > 0 && "kmax must be positive");
        blasfeo_dcolpei(kmax, ipiv, sA);
    }

    static inline void blasfeo_drowsw_debug(int kmax,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dmat *sC, int ci, int cj)
    {
        fatrop_dbg_assert(kmax > 0 && "kmax must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && ci >= 0 && cj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(aj + kmax <= sA->m && cj + kmax <= sC->m && "Row indices plus kmax must not exceed matrix dimensions");
        blasfeo_drowsw(kmax, sA, ai, aj, sC, ci, cj);
    }

    static inline void blasfeo_dcolsw_debug(int kmax,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dmat *sC, int ci, int cj)
    {
        fatrop_dbg_assert(kmax > 0 && "kmax must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && ci >= 0 && cj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + kmax <= sA->n && ci + kmax <= sC->n && "Column indices plus kmax must not exceed matrix dimensions");
        blasfeo_dcolsw(kmax, sA, ai, aj, sC, ci, cj);
    }

    static inline void blasfeo_dgead_debug(int m, int n, double alpha,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dmat *sB, int bi, int bj)
    {
        fatrop_dbg_assert(m > 0 && n > 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n && bi + m <= sB->m && bj + n <= sB->n && "Submatrix must fit within matrix dimensions");
        blasfeo_dgead(m, n, alpha, sA, ai, aj, sB, bi, bj);
    }

    static inline void blasfeo_dgecp_debug(int m, int n,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dmat *sB, int bi, int bj)
    {
        fatrop_dbg_assert(m > 0 && n > 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n && bi + m <= sB->m && bj + n <= sB->n && "Submatrix must fit within matrix dimensions");
        blasfeo_dgecp(m, n, sA, ai, aj, sB, bi, bj);
    }

    static inline void blasfeo_dgesc_debug(int m, int n, double alpha,  blasfeo_dmat *sA, int ai, int aj)
    {
        fatrop_dbg_assert(m > 0 && n > 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n && "Submatrix must fit within matrix dimensions");
        blasfeo_dgesc(m, n, alpha, sA, ai, aj);
    }

    static inline void blasfeo_dtrsm_rltn_debug(int m, int n, double alpha,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dmat *sB, int bi, int bj,  blasfeo_dmat *sD, int di, int dj)
    {
        fatrop_dbg_assert(m > 0 && n > 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 && di >= 0 && dj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + m <= sA->n && bi + m <= sB->m && bj + n <= sB->n && di + m <= sD->m && dj + n <= sD->n && "Submatrices must fit within matrix dimensions");
        blasfeo_dtrsm_rltn(m, n, alpha, sA, ai, aj, sB, bi, bj, sD, di, dj);
    }

    static inline void blasfeo_dgemm_nt_debug(int m, int n, int k, double alpha,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dmat *sB, int bi, int bj, double beta,  blasfeo_dmat *sC, int ci, int cj,  blasfeo_dmat *sD, int di, int dj)
    {
        fatrop_dbg_assert(m > 0 && n > 0 && k > 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 && ci >= 0 && cj >= 0 && di >= 0 && dj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + k <= sA->n && bi + n <= sB->m && bj + k <= sB->n && ci + m <= sC->m && cj + n <= sC->n && di + m <= sD->m && dj + n <= sD->n && "Submatrices must fit within matrix dimensions");
        blasfeo_dgemm_nt(m, n, k, alpha, sA, ai, aj, sB, bi, bj, beta, sC, ci, cj, sD, di, dj);
    }

    static inline void blasfeo_dsyrk_ln_mn_debug(int m, int n, int k, double alpha,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dmat *sB, int bi, int bj, double beta,  blasfeo_dmat *sC, int ci, int cj,  blasfeo_dmat *sD, int di, int dj)
    {
        fatrop_dbg_assert(m > 0 && n > 0 && k > 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 && ci >= 0 && cj >= 0 && di >= 0 && dj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + k <= sA->n && bi + n <= sB->m && bj + k <= sB->n && ci + m <= sC->m && cj + n <= sC->n && di + m <= sD->m && dj + n <= sD->n && "Submatrices must fit within matrix dimensions");
        blasfeo_dsyrk_ln_mn(m, n, k, alpha, sA, ai, aj, sB, bi, bj, beta, sC, ci, cj, sD, di, dj);
    }

    static inline void blasfeo_dsyrk_ln_debug(int m, int k, double alpha,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dmat *sB, int bi, int bj, double beta,  blasfeo_dmat *sC, int ci, int cj,  blasfeo_dmat *sD, int di, int dj)
    {
        fatrop_dbg_assert(m > 0 && k > 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 && ci >= 0 && cj >= 0 && di >= 0 && dj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + k <= sA->n && bi + m <= sB->m && bj + k <= sB->n && ci + m <= sC->m && cj + m <= sC->n && di + m <= sD->m && dj + m <= sD->n && "Submatrices must fit within matrix dimensions");
        blasfeo_dsyrk_ln(m, k, alpha, sA, ai, aj, sB, bi, bj, beta, sC, ci, cj, sD, di, dj);
    }

    static inline void blasfeo_dgetr_debug(int m, int n,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dmat *sB, int bi, int bj)
    {
        fatrop_dbg_assert(m > 0 && n > 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n && bi + n <= sB->m && bj + m <= sB->n && "Submatrices must fit within matrix dimensions");
        blasfeo_dgetr(m, n, sA, ai, aj, sB, bi, bj);
    }

    static inline void blasfeo_dtrtr_l_debug(int m,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dmat *sB, int bi, int bj)
    {
        fatrop_dbg_assert(m > 0 && "Matrix dimension must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + m <= sA->n && bi + m <= sB->m && bj + m <= sB->n && "Submatrices must fit within matrix dimensions");
        blasfeo_dtrtr_l(m, sA, ai, aj, sB, bi, bj);
    }

    static inline void blasfeo_dpotrf_l_mn_debug(int m, int n,  blasfeo_dmat *sC, int ci, int cj,  blasfeo_dmat *sD, int di, int dj)
    {
        fatrop_dbg_assert(m > 0 && n > 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ci >= 0 && cj >= 0 && di >= 0 && dj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ci + m <= sC->m && cj + n <= sC->n && di + m <= sD->m && dj + n <= sD->n && "Submatrices must fit within matrix dimensions");
        blasfeo_dpotrf_l_mn(m, n, sC, ci, cj, sD, di, dj);
    }

    static inline void blasfeo_drowex_debug(int kmax, double alpha,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dvec *sx, int xi)
    {
        fatrop_dbg_assert(kmax > 0 && "kmax must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && xi >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai < sA->m && aj + kmax <= sA->n && xi + kmax <= sx->m && "Indices and kmax must fit within matrix and vector dimensions");
        blasfeo_drowex(kmax, alpha, sA, ai, aj, sx, xi);
    }

    static inline void blasfeo_drowin_debug(int kmax, double alpha, blasfeo_dvec *sx, int xi,  blasfeo_dmat *sA, int ai, int aj)
    {
        fatrop_dbg_assert(kmax > 0 && "kmax must be positive");
        fatrop_dbg_assert(xi >= 0 && ai >= 0 && aj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(xi + kmax <= sx->m && ai < sA->m && aj + kmax <= sA->n && "Indices and kmax must fit within vector and matrix dimensions");
        blasfeo_drowin(kmax, alpha, sx, xi, sA, ai, aj);
    }

    static inline void blasfeo_dcolin_debug(int kmax,  blasfeo_dvec *sx, int xi,  blasfeo_dmat *sA, int ai, int aj)
    {
        fatrop_dbg_assert(kmax > 0 && "kmax must be positive");
        fatrop_dbg_assert(xi >= 0 && ai >= 0 && aj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(xi + kmax <= sx->m && ai + kmax <= sA->m && aj < sA->n && "Indices and kmax must fit within vector and matrix dimensions");
        blasfeo_dcolin(kmax, sx, xi, sA, ai, aj);
    }

    static inline void blasfeo_dtrsv_ltn_debug(int m,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dvec *sx, int xi,  blasfeo_dvec *sz, int zi)
    {
        fatrop_dbg_assert(m > 0 && "Matrix dimension must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && xi >= 0 && zi >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + m <= sA->n && xi + m <= sx->m && zi + m <= sz->m && "Indices and m must fit within matrix and vector dimensions");
        blasfeo_dtrsv_ltn(m, sA, ai, aj, sx, xi, sz, zi);
    }

    static inline void blasfeo_dtrsv_lnn_debug(int m,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dvec *sx, int xi,  blasfeo_dvec *sz, int zi)
    {
        fatrop_dbg_assert(m > 0 && "Matrix dimension must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && xi >= 0 && zi >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + m <= sA->n && xi + m <= sx->m && zi + m <= sz->m && "Indices and m must fit within matrix and vector dimensions");
        blasfeo_dtrsv_lnn(m, sA, ai, aj, sx, xi, sz, zi);
    }

    static inline void blasfeo_dtrsv_utn_debug(int m,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dvec *sx, int xi,  blasfeo_dvec *sz, int zi)
    {
        fatrop_dbg_assert(m > 0 && "Matrix dimension must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && xi >= 0 && zi >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + m <= sA->n && xi + m <= sx->m && zi + m <= sz->m && "Indices and m must fit within matrix and vector dimensions");
        blasfeo_dtrsv_utn(m, sA, ai, aj, sx, xi, sz, zi);
    }

    static inline void blasfeo_dgemv_t_debug(int m, int n, double alpha,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dvec *sx, int xi, double beta,  blasfeo_dvec *sy, int yi,  blasfeo_dvec *sz, int zi)
    {
        fatrop_dbg_assert(m > 0 && n > 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && xi >= 0 && yi >= 0 && zi >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n && xi + n <= sx->m && yi + m <= sy->m && zi + m <= sz->m && "Indices and dimensions must fit within matrix and vector dimensions");
        blasfeo_dgemv_t(m, n, alpha, sA, ai, aj, sx, xi, beta, sy, yi, sz, zi);
    }

    static inline void blasfeo_dgemv_n_debug(int m, int n, double alpha,  blasfeo_dmat *sA, int ai, int aj,  blasfeo_dvec *sx, int xi, double beta,  blasfeo_dvec *sy, int yi,  blasfeo_dvec *sz, int zi)
    {
        fatrop_dbg_assert(m > 0 && n > 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && xi >= 0 && yi >= 0 && zi >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n && xi + n <= sx->m && yi + m <= sy->m && zi + m <= sz->m && "Indices and dimensions must fit within matrix and vector dimensions");
        blasfeo_dgemv_n(m, n, alpha, sA, ai, aj, sx, xi, beta, sy, yi, sz, zi);
    }

    static inline void blasfeo_pack_dmat_debug(int m, int n, double *A, int lda,  blasfeo_dmat *sA, int ai, int aj)
    {
        fatrop_dbg_assert(m > 0 && n > 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(lda >= m && "lda must be greater than or equal to m");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n && "Indices and dimensions must fit within matrix dimensions");
        blasfeo_pack_dmat(m, n, A, lda, sA, ai, aj);
    }

    static inline void blasfeo_unpack_dvec_debug(int m,  blasfeo_dvec *sx, int xi, double *x, int incx)
    {
        fatrop_dbg_assert(m > 0 && "Vector dimension must be positive");
        fatrop_dbg_assert(xi >= 0 && "Index must be non-negative");
        fatrop_dbg_assert(xi + m <= sx->m && "Index and dimension must fit within vector dimension");
        fatrop_dbg_assert(incx != 0 && "incx must not be zero");
        blasfeo_unpack_dvec(m, sx, xi, x, incx);
    }

    static inline void blasfeo_pack_dvec_debug(int m, double *x, int incx,  blasfeo_dvec *sx, int xi)
    {
        fatrop_dbg_assert(m > 0 && "Vector dimension must be positive");
        fatrop_dbg_assert(xi >= 0 && "Index must be non-negative");
        fatrop_dbg_assert(xi + m <= sx->m && "Index and dimension must fit within vector dimension");
        fatrop_dbg_assert(incx != 0 && "incx must not be zero");
        blasfeo_pack_dvec(m, x, incx, sx, xi);
    }

    static inline double blasfeo_ddot_debug(int m,  blasfeo_dvec *sx, int xi,  blasfeo_dvec *sy, int yi)
    {
        fatrop_dbg_assert(m > 0 && "Vector dimension must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(xi + m <= sx->m && yi + m <= sy->m && "Indices and dimension must fit within vector dimensions");
        return blasfeo_ddot(m, sx, xi, sy, yi);
    }

    static inline void blasfeo_dgese_debug(int m, int n, double alpha, blasfeo_dmat *sA, int ai, int aj)
    {
        fatrop_dbg_assert(m > 0 && n > 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n && "Submatrix must fit within matrix dimensions");
        blasfeo_dgese(m, n, alpha, sA, ai, aj);
    }

    static inline void blasfeo_ddiare_debug(int kmax, double alpha, blasfeo_dmat *sA, int ai, int aj)
    {
        fatrop_dbg_assert(kmax > 0 && "kmax must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + kmax <= sA->m && aj + kmax <= sA->n && "Submatrix must fit within matrix dimensions");
        blasfeo_ddiare(kmax, alpha, sA, ai, aj);
    }

    static inline void blasfeo_dcolsc_debug(int kmax, double alpha, blasfeo_dmat *sA, int ai, int aj)
    {
        fatrop_dbg_assert(kmax > 0 && "kmax must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + kmax <= sA->m && aj < sA->n && "Column must fit within matrix dimensions");
        blasfeo_dcolsc(kmax, alpha, sA, ai, aj);
    }

    static inline void blasfeo_dvecmulacc_debug(int m, blasfeo_dvec *sx, int xi, blasfeo_dvec *sy, int yi, blasfeo_dvec *sz, int zi)
    {
        fatrop_dbg_assert(m > 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && zi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= sx->m && yi + m <= sy->m && zi + m <= sz->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        blasfeo_dvecmulacc(m, sx, xi, sy, yi, sz, zi);
    }

    static inline void blasfeo_dger_debug(int m, int n, double alpha, blasfeo_dvec *sx, int xi, blasfeo_dvec *sy, int yi, blasfeo_dmat *sC, int ci, int cj, blasfeo_dmat *sD, int di, int dj)
    {
        fatrop_dbg_assert(m > 0 && n > 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && ci >= 0 && cj >= 0 && di >= 0 && dj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(xi + m <= sx->m && yi + n <= sy->m && "Vector indices plus size must not exceed vector dimensions");
        fatrop_dbg_assert(ci + m <= sC->m && cj + n <= sC->n && di + m <= sD->m && dj + n <= sD->n && "Submatrices must fit within matrix dimensions");
        blasfeo_dger(m, n, alpha, sx, xi, sy, yi, sC, ci, cj, sD, di, dj);
    }

}

#endif // __fatrop_linear_algebra_blasfeo_wrapper_hpp__
