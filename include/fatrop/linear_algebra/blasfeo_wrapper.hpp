//
// Copyright (c) Lander Vanroye, KU Leuven
//
#ifndef __fatrop_linear_algebra_blasfeo_wrapper_hpp__
#define __fatrop_linear_algebra_blasfeo_wrapper_hpp__

/**
 * @file blasfeo_wrapper.hpp
 * @brief The functions in this file wrap the blasfeo functions. The wrapping adds two things not
 * provided by blasfeo by default (1) bound checking (in debug mode), and (2) correct constness, we
 * use const_casts to this end*/

#include "fatrop/common/exception.hpp"
extern "C"
{
#include <blasfeo.h>
}
#include "fatrop/context/context.hpp"

namespace fatrop
{
    static inline Scalar &blasfeo_matel_wrap(MAT *A, const int ai, const int aj)
    {
        fatrop_dbg_assert(ai < A->m && aj < A->n);
        return MATEL(A, ai, aj);
    }
    static inline Scalar &blasfeo_vecel_wrap(VEC *A, const int ai)
    {
        fatrop_dbg_assert(ai < A->m);
        return VECEL(A, ai);
    }

    static inline Scalar blasfeo_matel_wrap(const MAT *A, const int ai, const int aj)
    {
        fatrop_dbg_assert(ai < A->m && aj < A->n);
        return MATEL(A, ai, aj);
    }
    static inline Scalar blasfeo_vecel_wrap(const VEC *A, const int ai)
    {
        fatrop_dbg_assert(ai < A->m);
        return VECEL(A, ai);
    }
    static inline void blasfeo_axpby_wrap(int m, Scalar alpha, const VEC *x, int xi, Scalar beta,
                                          const VEC *y, int yi, VEC *z, int zi)
    {
        fatrop_dbg_assert(m >= 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && zi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m && yi + m <= y->m && zi + m <= z->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        AXPBY(m, alpha, const_cast<VEC *>(x), xi, beta, const_cast<VEC *>(y), yi, z, zi);
    }

    static inline void blasfeo_axpy_wrap(int m, Scalar alpha, const VEC *x, int xi, const VEC *y,
                                         int yi, VEC *z, int zi)
    {
        fatrop_dbg_assert(m >= 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && zi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m && yi + m <= y->m && zi + m <= z->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        AXPY(m, alpha, const_cast<VEC *>(x), xi, const_cast<VEC *>(y), yi, z, zi);
    }

    static inline void blasfeo_vecse_wrap(int m, Scalar alpha, VEC *x, int xi)
    {
        fatrop_dbg_assert(m >= 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && "Vector index must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m &&
                          "Vector index plus size must not exceed vector dimension");
        VECSE(m, alpha, x, xi);
    }

    static inline void blasfeo_vecsc_wrap(int m, Scalar alpha, VEC *x, int xi)
    {
        fatrop_dbg_assert(m >= 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && "Vector index must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m &&
                          "Vector index plus size must not exceed vector dimension");
        VECSC(m, alpha, x, xi);
    }

    static inline void blasfeo_veccpsc_wrap(int m, Scalar alpha, const VEC *x, int xi, VEC *y,
                                            int yi)
    {
        fatrop_dbg_assert(m >= 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m && yi + m <= y->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        VECCPSC(m, alpha, const_cast<VEC *>(x), xi, y, yi);
    }

    static inline void blasfeo_veccp_wrap(int m, const VEC *x, int xi, VEC *y, int yi)
    {
        fatrop_dbg_assert(m >= 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m && yi + m <= y->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        VECCP(m, const_cast<VEC *>(x), xi, y, yi);
    }

    static inline void blasfeo_vecmul_wrap(int m, const VEC *x, int xi, const VEC *y, int yi,
                                           VEC *z, int zi)
    {
        fatrop_dbg_assert(m >= 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && zi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= x->m && yi + m <= y->m && zi + m <= z->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        VECMUL(m, const_cast<VEC *>(x), xi, const_cast<VEC *>(y), yi, z, zi);
    }

    static inline void blasfeo_rowpe_wrap(int kmax, int *ipiv, MAT *sA)
    {
        fatrop_dbg_assert(kmax >= 0 && "kmax must be positive");
        ROWPE(kmax, ipiv, sA);
    }

    static inline void blasfeo_vecpe_wrap(int kmax, int *ipiv, VEC *sx, int xi)
    {
        fatrop_dbg_assert(kmax >= 0 && "kmax must be positive");
        fatrop_dbg_assert(xi >= 0 && "xi must be non-negative");
        fatrop_dbg_assert(xi + kmax <= sx->m && "xi + kmax must not exceed vector dimension");
        VECPE(kmax, ipiv, sx, xi);
    }

    static inline void blasfeo_vecpei_wrap(int kmax, int *ipiv, VEC *sx, int xi)
    {
        fatrop_dbg_assert(kmax >= 0 && "kmax must be positive");
        fatrop_dbg_assert(xi >= 0 && "xi must be non-negative");
        fatrop_dbg_assert(xi + kmax <= sx->m && "xi + kmax must not exceed vector dimension");
        VECPEI(kmax, ipiv, sx, xi);
    }

    static inline void blasfeo_rowpei_wrap(int kmax, int *ipiv, MAT *sA)
    {
        fatrop_dbg_assert(kmax >= 0 && "kmax must be positive");
        ROWPEI(kmax, ipiv, sA);
    }

    static inline void blasfeo_colpe_wrap(int kmax, int *ipiv, MAT *sA)
    {
        fatrop_dbg_assert(kmax >= 0 && "kmax must be positive");
        COLPE(kmax, ipiv, sA);
    }

    static inline void blasfeo_colpei_wrap(int kmax, int *ipiv, MAT *sA)
    {
        fatrop_dbg_assert(kmax >= 0 && "kmax must be positive");
        COLPEI(kmax, ipiv, sA);
    }

    static inline void blasfeo_rowsw_wrap(int kmax, MAT *sA, int ai, int aj, MAT *sC, int ci,
                                          int cj)
    {
        fatrop_dbg_assert(kmax >= 0 && "kmax must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && ci >= 0 && cj >= 0 &&
                          "Indices must be non-negative");
        fatrop_dbg_assert(aj + kmax <= sA->m && cj + kmax <= sC->m &&
                          "Row indices plus kmax must not exceed matrix dimensions");
        ROWSW(kmax, sA, ai, aj, sC, ci, cj);
    }

    static inline void blasfeo_colsw_wrap(int kmax, MAT *sA, int ai, int aj, MAT *sC, int ci,
                                          int cj)
    {
        fatrop_dbg_assert(kmax >= 0 && "kmax must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && ci >= 0 && cj >= 0 &&
                          "Indices must be non-negative");
        fatrop_dbg_assert(ai + kmax <= sA->n && ci + kmax <= sC->n &&
                          "Column indices plus kmax must not exceed matrix dimensions");
        COLSW(kmax, sA, ai, aj, sC, ci, cj);
    }

    static inline void blasfeo_gead_wrap(int m, int n, Scalar alpha, const MAT *sA, int ai, int aj,
                                         MAT *sB, int bi, int bj)
    {
        fatrop_dbg_assert(m >= 0 && n >= 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 &&
                          "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n && bi + m <= sB->m &&
                          bj + n <= sB->n && "Submatrix must fit within matrix dimensions");
        GEAD(m, n, alpha, const_cast<MAT *>(sA), ai, aj, sB, bi, bj);
    }

    static inline void blasfeo_gecp_wrap(int m, int n, const MAT *sA, int ai, int aj, MAT *sB,
                                         int bi, int bj)
    {
        fatrop_dbg_assert(m >= 0 && n >= 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 &&
                          "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n && bi + m <= sB->m &&
                          bj + n <= sB->n && "Submatrix must fit within matrix dimensions");
        GECP(m, n, const_cast<MAT *>(sA), ai, aj, sB, bi, bj);
    }

    static inline void blasfeo_gesc_wrap(int m, int n, Scalar alpha, MAT *sA, int ai, int aj)
    {
        fatrop_dbg_assert(m >= 0 && n >= 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n &&
                          "Submatrix must fit within matrix dimensions");
        GESC(m, n, alpha, sA, ai, aj);
    }

    static inline void blasfeo_trsm_rltn_wrap(int m, int n, Scalar alpha, const MAT *sA, int ai,
                                              int aj, const MAT *sB, int bi, int bj, MAT *sD,
                                              int di, int dj)
    {
        // D <= alpha * B * A^{-T} , with A lower triangular employing explicit inverse of diagonal
        // B and D mxn, A nxn
        fatrop_dbg_assert(m >= 0 && n >= 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 && di >= 0 && dj >= 0 &&
                          "Indices must be non-negative");
        fatrop_dbg_assert(ai + n <= sA->m && aj + n <= sA->n && bi + m <= sB->m &&
                          bj + n <= sB->n && di + m <= sD->m && dj + n <= sD->n &&
                          "Submatrices must fit within matrix dimensions");
        TRSM_RLTN(m, n, alpha, const_cast<MAT *>(sA), ai, aj, const_cast<MAT *>(sB), bi, bj, sD, di,
                  dj);
    }

    static inline void blasfeo_trsm_rlnn_wrap(int m, int n, Scalar alpha, const MAT *sA, int ai,
                                              int aj, const MAT *sB, int bi, int bj, MAT *sD,
                                              int di, int dj)
    {
        fatrop_dbg_assert(m >= 0 && n >= 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 && di >= 0 && dj >= 0 &&
                          "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + m <= sA->n && bi + m <= sB->m &&
                          bj + n <= sB->n && di + m <= sD->m && dj + n <= sD->n &&
                          "Submatrices must fit within matrix dimensions");
        TRSM_RLNN(m, n, alpha, const_cast<MAT *>(sA), ai, aj, const_cast<MAT *>(sB), bi, bj, sD, di,
                  dj);
    }

    static inline void blasfeo_gemm_nt_wrap(int m, int n, int k, Scalar alpha, const MAT *sA,
                                            int ai, int aj, const MAT *sB, int bi, int bj,
                                            Scalar beta, const MAT *sC, int ci, int cj, MAT *sD,
                                            int di, int dj)
    {
        fatrop_dbg_assert(m >= 0 && n >= 0 && k >= 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 && ci >= 0 && cj >= 0 &&
                          di >= 0 && dj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + k <= sA->n && bi + n <= sB->m &&
                          bj + k <= sB->n && ci + m <= sC->m && cj + n <= sC->n &&
                          di + m <= sD->m && dj + n <= sD->n &&
                          "Submatrices must fit within matrix dimensions");
        GEMM_NT(m, n, k, alpha, const_cast<MAT *>(sA), ai, aj, const_cast<MAT *>(sB), bi, bj, beta,
                const_cast<MAT *>(sC), ci, cj, sD, di, dj);
    }

    static inline void blasfeo_syrk_ln_mn_wrap(int m, int n, int k, Scalar alpha, const MAT *sA,
                                               int ai, int aj, const MAT *sB, int bi, int bj,
                                               Scalar beta, const MAT *sC, int ci, int cj, MAT *sD,
                                               int di, int dj)
    {
        fatrop_dbg_assert(m >= 0 && n >= 0 && k >= 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 && ci >= 0 && cj >= 0 &&
                          di >= 0 && dj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + k <= sA->n && bi + n <= sB->m &&
                          bj + k <= sB->n && ci + m <= sC->m && cj + n <= sC->n &&
                          di + m <= sD->m && dj + n <= sD->n &&
                          "Submatrices must fit within matrix dimensions");
        SYRK_LN_MN(m, n, k, alpha, const_cast<MAT *>(sA), ai, aj, const_cast<MAT *>(sB), bi, bj,
                   beta, const_cast<MAT *>(sC), ci, cj, sD, di, dj);
    }

    static inline void blasfeo_syrk_ln_wrap(int m, int k, Scalar alpha, const MAT *sA, int ai,
                                            int aj, const MAT *sB, int bi, int bj, Scalar beta,
                                            const MAT *sC, int ci, int cj, MAT *sD, int di, int dj)
    {
        fatrop_dbg_assert(m >= 0 && k >= 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 && ci >= 0 && cj >= 0 &&
                          di >= 0 && dj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + k <= sA->n && bi + m <= sB->m &&
                          bj + k <= sB->n && ci + m <= sC->m && cj + m <= sC->n &&
                          di + m <= sD->m && dj + m <= sD->n &&
                          "Submatrices must fit within matrix dimensions");
        SYRK_LN(m, k, alpha, const_cast<MAT *>(sA), ai, aj, const_cast<MAT *>(sB), bi, bj, beta,
                const_cast<MAT *>(sC), ci, cj, sD, di, dj);
    }

    static inline void blasfeo_getr_wrap(int m, int n, const MAT *sA, int ai, int aj, MAT *sB,
                                         int bi, int bj)
    {
        fatrop_dbg_assert(m >= 0 && n >= 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 &&
                          "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n && bi + n <= sB->m &&
                          bj + m <= sB->n && "Submatrices must fit within matrix dimensions");
        GETR(m, n, const_cast<MAT *>(sA), ai, aj, sB, bi, bj);
    }

    static inline void blasfeo_trtr_l_wrap(int m, const MAT *sA, int ai, int aj, MAT *sB, int bi,
                                           int bj)
    {
        fatrop_dbg_assert(m >= 0 && "Matrix dimension must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && bi >= 0 && bj >= 0 &&
                          "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + m <= sA->n && bi + m <= sB->m &&
                          bj + m <= sB->n && "Submatrices must fit within matrix dimensions");
        TRTR_L(m, const_cast<MAT *>(sA), ai, aj, sB, bi, bj);
    }

    static inline void blasfeo_potrf_l_mn_wrap(int m, int n, MAT *sC, int ci, int cj, MAT *sD,
                                               int di, int dj)
    {
        fatrop_dbg_assert(m >= 0 && n >= 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ci >= 0 && cj >= 0 && di >= 0 && dj >= 0 &&
                          "Indices must be non-negative");
        fatrop_dbg_assert(ci + m <= sC->m && cj + n <= sC->n && di + m <= sD->m &&
                          dj + n <= sD->n && "Submatrices must fit within matrix dimensions");
        POTRF_L_MN(m, n, sC, ci, cj, sD, di, dj);
    }

    static inline void blasfeo_rowex_wrap(int kmax, Scalar alpha, const MAT *sA, int ai, int aj,
                                          VEC *sx, int xi)
    {
        fatrop_dbg_assert(kmax >= 0 && "kmax must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && xi >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai < sA->m && aj + kmax <= sA->n && xi + kmax <= sx->m &&
                          "Indices and kmax must fit within matrix and vector dimensions");
        ROWEX(kmax, alpha, const_cast<MAT *>(sA), ai, aj, sx, xi);
    }

    static inline void blasfeo_rowin_wrap(int kmax, Scalar alpha, const VEC *sx, int xi, MAT *sA,
                                          int ai, int aj)
    {
        fatrop_dbg_assert(kmax >= 0 && "kmax must be positive");
        fatrop_dbg_assert(xi >= 0 && ai >= 0 && aj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(xi + kmax <= sx->m && ai < sA->m && aj + kmax <= sA->n &&
                          "Indices and kmax must fit within vector and matrix dimensions");
        ROWIN(kmax, alpha, const_cast<VEC *>(sx), xi, sA, ai, aj);
    }

    static inline void blasfeo_colin_wrap(int kmax, const VEC *sx, int xi, MAT *sA, int ai, int aj)
    {
        fatrop_dbg_assert(kmax >= 0 && "kmax must be positive");
        fatrop_dbg_assert(xi >= 0 && ai >= 0 && aj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(xi + kmax <= sx->m && ai + kmax <= sA->m && aj < sA->n &&
                          "Indices and kmax must fit within vector and matrix dimensions");
        COLIN(kmax, const_cast<VEC *>(sx), xi, sA, ai, aj);
    }

    static inline void blasfeo_trsv_ltn_wrap(int m, const MAT *sA, int ai, int aj, const VEC *sx,
                                             int xi, VEC *sz, int zi)
    {
        fatrop_dbg_assert(m >= 0 && "Matrix dimension must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && xi >= 0 && zi >= 0 &&
                          "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + m <= sA->n && xi + m <= sx->m &&
                          zi + m <= sz->m &&
                          "Indices and m must fit within matrix and vector dimensions");
        TRSV_LTN(m, const_cast<MAT *>(sA), ai, aj, const_cast<VEC *>(sx), xi, sz, zi);
    }

    static inline void blasfeo_trsv_lnn_wrap(int m, const MAT *sA, int ai, int aj, const VEC *sx,
                                             int xi, VEC *sz, int zi)
    {
        fatrop_dbg_assert(m >= 0 && "Matrix dimension must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && xi >= 0 && zi >= 0 &&
                          "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + m <= sA->n && xi + m <= sx->m &&
                          zi + m <= sz->m &&
                          "Indices and m must fit within matrix and vector dimensions");
        TRSV_LNN(m, const_cast<MAT *>(sA), ai, aj, const_cast<VEC *>(sx), xi, sz, zi);
    }

    static inline void blasfeo_trsv_utn_wrap(int m, const MAT *sA, int ai, int aj, const VEC *sx,
                                             int xi, VEC *sz, int zi)
    {
        fatrop_dbg_assert(m >= 0 && "Matrix dimension must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && xi >= 0 && zi >= 0 &&
                          "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + m <= sA->n && xi + m <= sx->m &&
                          zi + m <= sz->m &&
                          "Indices and m must fit within matrix and vector dimensions");
        TRSV_UTN(m, const_cast<MAT *>(sA), ai, aj, const_cast<VEC *>(sx), xi, sz, zi);
    }

    static inline void blasfeo_gemv_t_wrap(int m, int n, Scalar alpha, const MAT *sA, int ai,
                                           int aj, const VEC *sx, int xi, Scalar beta,
                                           const VEC *sy, int yi, VEC *sz, int zi)
    {
        fatrop_dbg_assert(m >= 0 && n >= 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && xi >= 0 && yi >= 0 && zi >= 0 &&
                          "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n && xi + m <= sx->m &&
                          yi + n <= sy->m && zi + n <= sz->m &&
                          "Indices and dimensions must fit within matrix and vector dimensions");
        GEMV_T(m, n, alpha, const_cast<MAT *>(sA), ai, aj, const_cast<VEC *>(sx), xi, beta,
               const_cast<VEC *>(sy), yi, sz, zi);
    }

    static inline void blasfeo_gemv_n_wrap(int m, int n, Scalar alpha, const MAT *sA, int ai,
                                           int aj, const VEC *sx, int xi, Scalar beta,
                                           const VEC *sy, int yi, VEC *sz, int zi)
    {
        fatrop_dbg_assert(m >= 0 && n >= 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && xi >= 0 && yi >= 0 && zi >= 0 &&
                          "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n && xi + n <= sx->m &&
                          yi + m <= sy->m && zi + m <= sz->m &&
                          "Indices and dimensions must fit within matrix and vector dimensions");
        GEMV_N(m, n, alpha, const_cast<MAT *>(sA), ai, aj, const_cast<VEC *>(sx), xi, beta,
               const_cast<VEC *>(sy), yi, sz, zi);
    }

    static inline void blasfeo_pack_mat_wrap(int m, int n, Scalar *A, int lda, MAT *sA, int ai,
                                             int aj)
    {
        fatrop_dbg_assert(m >= 0 && n >= 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(lda >= m && "lda must be greater than or equal to m");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n &&
                          "Indices and dimensions must fit within matrix dimensions");
        PACK_MAT(m, n, A, lda, sA, ai, aj);
    }

    static inline void blasfeo_unpack_vec_wrap(int m, const VEC *sx, int xi, Scalar *x, int incx)
    {
        fatrop_dbg_assert(m >= 0 && "Vector dimension must be positive");
        fatrop_dbg_assert(xi >= 0 && "Index must be non-negative");
        fatrop_dbg_assert(xi + m <= sx->m &&
                          "Index and dimension must fit within vector dimension");
        fatrop_dbg_assert(incx != 0 && "incx must not be zero");
        UNPACK_VEC(m, const_cast<VEC *>(sx), xi, x, incx);
    }

    static inline void blasfeo_pack_vec_wrap(int m, Scalar *x, int incx, VEC *sx, int xi)
    {
        fatrop_dbg_assert(m >= 0 && "Vector dimension must be positive");
        fatrop_dbg_assert(xi >= 0 && "Index must be non-negative");
        fatrop_dbg_assert(xi + m <= sx->m &&
                          "Index and dimension must fit within vector dimension");
        fatrop_dbg_assert(incx != 0 && "incx must not be zero");
        PACK_VEC(m, x, incx, sx, xi);
    }

    static inline Scalar blasfeo_dot_wrap(int m, const VEC *sx, int xi, const VEC *sy, int yi)
    {
        fatrop_dbg_assert(m >= 0 && "Vector dimension must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(xi + m <= sx->m && yi + m <= sy->m &&
                          "Indices and dimension must fit within vector dimensions");
        return DOT(m, const_cast<VEC *>(sx), xi, const_cast<VEC *>(sy), yi);
    }

    static inline void blasfeo_gese_wrap(int m, int n, Scalar alpha, MAT *sA, int ai, int aj)
    {
        fatrop_dbg_assert(m >= 0 && n >= 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + m <= sA->m && aj + n <= sA->n &&
                          "Submatrix must fit within matrix dimensions");
        GESE(m, n, alpha, sA, ai, aj);
    }

    static inline void blasfeo_diare_wrap(int kmax, Scalar alpha, MAT *sA, int ai, int aj)
    {
        fatrop_dbg_assert(kmax >= 0 && "kmax must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + kmax <= sA->m && aj + kmax <= sA->n &&
                          "Submatrix must fit within matrix dimensions");
        DIARE(kmax, alpha, sA, ai, aj);
    }

    static inline void blasfeo_diaad_wrap(int kmax, Scalar alpha, const VEC *sx, int xi, MAT *sA,
                                          int ai, int aj)
    {
        fatrop_dbg_assert(kmax >= 0 && "kmax must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && xi >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(xi + kmax <= sx -> m && "Subvector must fit within vector dimensions");
        fatrop_dbg_assert(ai + kmax <= sA->m && aj + kmax <= sA->n &&
                          "Submatrix must fit within matrix dimensions");
        DIAAD(kmax, alpha, const_cast<VEC *>(sx), xi, sA, ai, aj);
    }

    static inline void blasfeo_colsc_wrap(int kmax, Scalar alpha, MAT *sA, int ai, int aj)
    {
        fatrop_dbg_assert(kmax >= 0 && "kmax must be positive");
        fatrop_dbg_assert(ai >= 0 && aj >= 0 && "Indices must be non-negative");
        fatrop_dbg_assert(ai + kmax <= sA->m && aj < sA->n &&
                          "Column must fit within matrix dimensions");
        COLSC(kmax, alpha, sA, ai, aj);
    }

    static inline void blasfeo_vecmulacc_wrap(int m, const VEC *sx, int xi, const VEC *sy, int yi,
                                              VEC *sz, int zi)
    {
        fatrop_dbg_assert(m >= 0 && "Vector size must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && zi >= 0 && "Vector indices must be non-negative");
        fatrop_dbg_assert(xi + m <= sx->m && yi + m <= sy->m && zi + m <= sz->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        VECMULACC(m, const_cast<VEC *>(sx), xi, const_cast<VEC *>(sy), yi, sz, zi);
    }

    static inline void blasfeo_ger_wrap(int m, int n, Scalar alpha, const VEC *sx, int xi,
                                        const VEC *sy, int yi, const MAT *sC, int ci, int cj,
                                        MAT *sD, int di, int dj)
    {
        fatrop_dbg_assert(m >= 0 && n >= 0 && "Matrix dimensions must be positive");
        fatrop_dbg_assert(xi >= 0 && yi >= 0 && ci >= 0 && cj >= 0 && di >= 0 && dj >= 0 &&
                          "Indices must be non-negative");
        fatrop_dbg_assert(xi + m <= sx->m && yi + n <= sy->m &&
                          "Vector indices plus size must not exceed vector dimensions");
        fatrop_dbg_assert(ci + m <= sC->m && cj + n <= sC->n && di + m <= sD->m &&
                          dj + n <= sD->n && "Submatrices must fit within matrix dimensions");
        GER(m, n, alpha, const_cast<VEC *>(sx), xi, const_cast<VEC *>(sy), yi,
            const_cast<MAT *>(sC), ci, cj, sD, di, dj);
    }

}

#endif // __fatrop_linear_algebra_blasfeo_wrapper_hpp__
