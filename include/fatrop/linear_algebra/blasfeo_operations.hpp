//
// Copyright (c) Lander Vanroye, KU Leuven
//
#ifndef __fatrop_linear_algebra_blasfeo_operations_hpp__
#define __fatrop_linear_algebra_blasfeo_operations_hpp__
/**
 * @file blasfeo_operations.hpp
 * @brief The functions in this file makes the blasfeo functions able to work with the fatrop types
 *for representing vectors and matrices.
 **/

#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/blasfeo_wrapper.hpp"
#include "fatrop/linear_algebra/lu_factorization.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/linear_algebra/vector.hpp"

namespace fatrop
{

    static inline void axpby(int m, Scalar alpha, const VecRealView &x, int xi, Scalar beta,
                             const VecRealView &y, int yi, VecRealView &z, int zi)
    {
        blasfeo_axpby_wrap(m, alpha, &x.vec(), x.ai() + xi, beta, &y.vec(), y.ai() + yi, &z.vec(),
                           z.ai() + zi);
    }

    static inline void axpy(int m, Scalar alpha, const VecRealView &x, int xi, const VecRealView &y,
                            int yi, VecRealView &z, int zi)
    {
        blasfeo_axpy_wrap(m, alpha, &x.vec(), x.ai() + xi, &y.vec(), y.ai() + yi, &z.vec(),
                          z.ai() + zi);
    }

    static inline void vecse(int m, Scalar alpha, VecRealView &x, int xi)
    {
        blasfeo_vecse_wrap(m, alpha, &x.vec(), x.ai() + xi);
    }

    static inline void vecsc(int m, Scalar alpha, VecRealView &x, int xi)
    {
        blasfeo_vecsc_wrap(m, alpha, &x.vec(), x.ai() + xi);
    }

    static inline void veccpsc(int m, Scalar alpha, const VecRealView &x, int xi, VecRealView &y,
                               int yi)
    {
        blasfeo_veccpsc_wrap(m, alpha, &x.vec(), x.ai() + xi, &y.vec(), y.ai() + yi);
    }

    static inline void veccp(int m, const VecRealView &x, int xi, VecRealView &y, int yi)
    {
        blasfeo_veccp_wrap(m, &x.vec(), x.ai() + xi, &y.vec(), y.ai() + yi);
    }

    static inline void vecmul(int m, const VecRealView &x, int xi, const VecRealView &y, int yi,
                              VecRealView &z, int zi)
    {
        blasfeo_vecmul_wrap(m, &x.vec(), x.ai() + xi, &y.vec(), y.ai() + yi, &z.vec(), z.ai() + zi);
    }

    static inline void rowpe(int kmax, int *ipiv, MatRealView &A)
    {
        blasfeo_rowpe_wrap(kmax, ipiv, &A.mat());
    }

    static inline void vecpe(int kmax, int *ipiv, VecRealView &x, int xi)
    {
        blasfeo_vecpe_wrap(kmax, ipiv, &x.vec(), x.ai() + xi);
    }

    static inline void vecpei(int kmax, int *ipiv, VecRealView &x, int xi)
    {
        blasfeo_vecpei_wrap(kmax, ipiv, &x.vec(), x.ai() + xi);
    }

    static inline void rowpei(int kmax, int *ipiv, MatRealView &A)
    {
        blasfeo_rowpei_wrap(kmax, ipiv, &A.mat());
    }

    static inline void colpe(int kmax, int *ipiv, MatRealView &A)
    {
        blasfeo_colpe_wrap(kmax, ipiv, &A.mat());
    }

    static inline void colpei(int kmax, int *ipiv, MatRealView &A)
    {
        blasfeo_colpei_wrap(kmax, ipiv, &A.mat());
    }

    static inline void rowsw(int kmax, MatRealView &A, int ai, int aj, MatRealView &C, int ci,
                             int cj)
    {
        blasfeo_rowsw_wrap(kmax, &A.mat(), A.ai() + ai, A.aj() + aj, &C.mat(), C.ai() + ci,
                           C.aj() + cj);
    }

    static inline void colsw(int kmax, MatRealView &A, int ai, int aj, MatRealView &C, int ci,
                             int cj)
    {
        blasfeo_colsw_wrap(kmax, &A.mat(), A.ai() + ai, A.aj() + aj, &C.mat(), C.ai() + ci,
                           C.aj() + cj);
    }

    static inline void gead(int m, int n, Scalar alpha, const MatRealView &A, int ai, int aj,
                            MatRealView &B, int bi, int bj)
    {
        blasfeo_gead_wrap(m, n, alpha, &A.mat(), A.ai() + ai, A.aj() + aj, &B.mat(), B.ai() + bi,
                          B.aj() + bj);
    }

    static inline void gecp(int m, int n, const MatRealView &A, int ai, int aj, MatRealView &B,
                            int bi, int bj)
    {
        blasfeo_gecp_wrap(m, n, &A.mat(), A.ai() + ai, A.aj() + aj, &B.mat(), B.ai() + bi,
                          B.aj() + bj);
    }

    static inline void gesc(int m, int n, Scalar alpha, MatRealView &A, int ai, int aj)
    {
        blasfeo_gesc_wrap(m, n, alpha, &A.mat(), A.ai() + ai, A.aj() + aj);
    }

    static inline void trsm_rltn(int m, int n, Scalar alpha, const MatRealView &A, int ai, int aj,
                                 const MatRealView &B, int bi, int bj, MatRealView &D, int di,
                                 int dj)
    {
        blasfeo_trsm_rltn_wrap(m, n, alpha, &A.mat(), A.ai() + ai, A.aj() + aj, &B.mat(),
                               B.ai() + bi, B.aj() + bj, &D.mat(), D.ai() + di, D.aj() + dj);
    }

    static inline void trsm_rlnn(int m, int n, Scalar alpha, const MatRealView &A, int ai, int aj,
                                 const MatRealView &B, int bi, int bj, MatRealView &D, int di,
                                 int dj)
    {
        blasfeo_trsm_rlnn_wrap(m, n, alpha, &A.mat(), A.ai() + ai, A.aj() + aj, &B.mat(),
                               B.ai() + bi, B.aj() + bj, &D.mat(), D.ai() + di, D.aj() + dj);
    }

    static inline void gemm_nt(int m, int n, int k, Scalar alpha, const MatRealView &A, int ai,
                               int aj, const MatRealView &B, int bi, int bj, Scalar beta,
                               const MatRealView &C, int ci, int cj, MatRealView &D, int di, int dj)
    {
        blasfeo_gemm_nt_wrap(m, n, k, alpha, &A.mat(), A.ai() + ai, A.aj() + aj, &B.mat(),
                             B.ai() + bi, B.aj() + bj, beta, &C.mat(), C.ai() + ci, C.aj() + cj,
                             &D.mat(), D.ai() + di, D.aj() + dj);
    }

    static inline void syrk_ln_mn(int m, int n, int k, Scalar alpha, const MatRealView &A, int ai,
                                  int aj, const MatRealView &B, int bi, int bj, Scalar beta,
                                  const MatRealView &C, int ci, int cj, MatRealView &D, int di,
                                  int dj)
    {
        blasfeo_syrk_ln_mn_wrap(m, n, k, alpha, &A.mat(), A.ai() + ai, A.aj() + aj, &B.mat(),
                                B.ai() + bi, B.aj() + bj, beta, &C.mat(), C.ai() + ci, C.aj() + cj,
                                &D.mat(), D.ai() + di, D.aj() + dj);
    }

    static inline void syrk_ln(int m, int k, Scalar alpha, const MatRealView &A, int ai, int aj,
                               const MatRealView &B, int bi, int bj, Scalar beta,
                               const MatRealView &C, int ci, int cj, MatRealView &D, int di, int dj)
    {
        blasfeo_syrk_ln_wrap(m, k, alpha, &A.mat(), A.ai() + ai, A.aj() + aj, &B.mat(), B.ai() + bi,
                             B.aj() + bj, beta, &C.mat(), C.ai() + ci, C.aj() + cj, &D.mat(),
                             D.ai() + di, D.aj() + dj);
    }

    static inline void getr(int m, int n, const MatRealView &A, int ai, int aj, MatRealView &B,
                            int bi, int bj)
    {
        blasfeo_getr_wrap(m, n, &A.mat(), A.ai() + ai, A.aj() + aj, &B.mat(), B.ai() + bi,
                          B.aj() + bj);
    }

    static inline void trtr_l(int m, const MatRealView &A, int ai, int aj, MatRealView &B, int bi,
                              int bj)
    {
        blasfeo_trtr_l_wrap(m, &A.mat(), A.ai() + ai, A.aj() + aj, &B.mat(), B.ai() + bi,
                            B.aj() + bj);
    }

    static inline void potrf_l_mn(int m, int n, MatRealView &C, int ci, int cj, MatRealView &D,
                                  int di, int dj)
    {
        blasfeo_potrf_l_mn_wrap(m, n, &C.mat(), C.ai() + ci, C.aj() + cj, &D.mat(), D.ai() + di,
                                D.aj() + dj);
    }

    static inline void rowex(int kmax, Scalar alpha, const MatRealView &A, int ai, int aj,
                             VecRealView &x, int xi)
    {
        blasfeo_rowex_wrap(kmax, alpha, &A.mat(), A.ai() + ai, A.aj() + aj, &x.vec(), x.ai() + xi);
    }

    static inline void rowin(int kmax, Scalar alpha, const VecRealView &x, int xi, MatRealView &A,
                             int ai, int aj)
    {
        blasfeo_rowin_wrap(kmax, alpha, &x.vec(), x.ai() + xi, &A.mat(), A.ai() + ai, A.aj() + aj);
    }

    static inline void colin(int kmax, const VecRealView &x, int xi, MatRealView &A, int ai, int aj)
    {
        blasfeo_colin_wrap(kmax, &x.vec(), x.ai() + xi, &A.mat(), A.ai() + ai, A.aj() + aj);
    }

    static inline void trsv_ltn(int m, const MatRealView &A, int ai, int aj, const VecRealView &x,
                                int xi, VecRealView &z, int zi)
    {
        blasfeo_trsv_ltn_wrap(m, &A.mat(), A.ai() + ai, A.aj() + aj, &x.vec(), x.ai() + xi,
                              &z.vec(), z.ai() + zi);
    }

    static inline void trsv_lnn(int m, const MatRealView &A, int ai, int aj, const VecRealView &x,
                                int xi, VecRealView &z, int zi)
    {
        blasfeo_trsv_lnn_wrap(m, &A.mat(), A.ai() + ai, A.aj() + aj, &x.vec(), x.ai() + xi,
                              &z.vec(), z.ai() + zi);
    }

    static inline void trsv_utn(int m, const MatRealView &A, int ai, int aj, const VecRealView &x,
                                int xi, VecRealView &z, int zi)
    {
        blasfeo_trsv_utn_wrap(m, &A.mat(), A.ai() + ai, A.aj() + aj, &x.vec(), x.ai() + xi,
                              &z.vec(), z.ai() + zi);
    }

    static inline void gemv_t(int m, int n, Scalar alpha, const MatRealView &A, int ai, int aj,
                              const VecRealView &x, int xi, Scalar beta, const VecRealView &y,
                              int yi, VecRealView &z, int zi)
    {
        blasfeo_gemv_t_wrap(m, n, alpha, &A.mat(), A.ai() + ai, A.aj() + aj, &x.vec(), x.ai() + xi,
                            beta, &y.vec(), y.ai() + yi, &z.vec(), z.ai() + zi);
    }

    static inline void gemv_n(int m, int n, Scalar alpha, const MatRealView &A, int ai, int aj,
                              const VecRealView &x, int xi, Scalar beta, const VecRealView &y,
                              int yi, VecRealView &z, int zi)
    {
        blasfeo_gemv_n_wrap(m, n, alpha, &A.mat(), A.ai() + ai, A.aj() + aj, &x.vec(), x.ai() + xi,
                            beta, &y.vec(), y.ai() + yi, &z.vec(), z.ai() + zi);
    }

    static inline void pack_mat(int m, int n, Scalar *A, int lda, MatRealView &sA, int ai, int aj)
    {
        blasfeo_pack_mat_wrap(m, n, A, lda, &sA.mat(), sA.ai() + ai, sA.aj() + aj);
    }

    static inline void unpack_vec(int m, const VecRealView &sx, int xi, Scalar *x, int incx)
    {
        blasfeo_unpack_vec_wrap(m, &sx.vec(), sx.ai() + xi, x, incx);
    }

    static inline void pack_vec(int m, Scalar *x, int incx, VecRealView &sx, int xi)
    {
        blasfeo_pack_vec_wrap(m, x, incx, &sx.vec(), sx.ai() + xi);
    }

    static inline Scalar dot(int m, const VecRealView &sx, int xi, const VecRealView &sy, int yi)
    {
        return blasfeo_dot_wrap(m, &sx.vec(), sx.ai() + xi, &sy.vec(), sy.ai() + yi);
    }

    static inline void gese(int m, int n, Scalar alpha, MatRealView &sA, int ai, int aj)
    {
        blasfeo_gese_wrap(m, n, alpha, &sA.mat(), sA.ai() + ai, sA.aj() + aj);
    }

    static inline void diare(int kmax, Scalar alpha, MatRealView &sA, int ai, int aj)
    {
        blasfeo_diare_wrap(kmax, alpha, &sA.mat(), sA.ai() + ai, sA.aj() + aj);
    }

    static inline void diaad(int kmax, Scalar alpha, const VecRealView &sx, int xi, MatRealView &sA,
                             int ai, int aj)
    {
        blasfeo_diaad_wrap(kmax, alpha, &sx.vec(), sx.ai() + xi, &sA.mat(), sA.ai() + ai, sA.aj() + aj);
    }

    static inline void colsc(int kmax, Scalar alpha, MatRealView &sA, int ai, int aj)
    {
        blasfeo_colsc_wrap(kmax, alpha, &sA.mat(), sA.ai() + ai, sA.aj() + aj);
    }

    static inline void vecmulacc(int m, const VecRealView &sx, int xi, const VecRealView &sy,
                                 int yi, VecRealView &sz, int zi)
    {
        blasfeo_vecmulacc_wrap(m, &sx.vec(), sx.ai() + xi, &sy.vec(), sy.ai() + yi, &sz.vec(),
                               sz.ai() + zi);
    }

    static inline void ger(int m, int n, Scalar alpha, const VecRealView &sx, int xi,
                           const VecRealView &sy, int yi, const MatRealView &sC, int ci, int cj,
                           MatRealView &sD, int di, int dj)
    {
        blasfeo_ger_wrap(m, n, alpha, &sx.vec(), sx.ai() + xi, &sy.vec(), sy.ai() + yi, &sC.mat(),
                         sC.ai() + ci, sC.aj() + cj, &sD.mat(), sD.ai() + di, sD.aj() + dj);
    }

    static inline void gead_transposed(int m, int n, Scalar alpha, const MatRealView &A, int ai,
                                       int aj, MatRealView &B, int bi, int bj)
    {
        fatrop_gead_transposed(m, n, alpha, const_cast<MAT *>(&A.mat()), A.ai() + ai, A.aj() + aj,
                               &B.mat(), B.ai() + bi, B.aj() + bj);
    }

    static inline void trsv_unu(int m, int n, const MatRealView &A, int ai, int aj,
                                const VecRealView &x, int xi, VecRealView &z, int zi)
    {
        fatrop_trsv_unu(m, n, const_cast<MAT *>(&A.mat()), A.ai() + ai, A.aj() + aj,
                        const_cast<VEC *>(&x.vec()), x.ai() + xi, &z.vec(), z.ai() + zi);
    }

    static inline void trsv_utu(int m, const MatRealView &A, int ai, int aj, const VecRealView &x,
                                int xi, VecRealView &z, int zi)
    {
        fatrop_trsv_utu(m, const_cast<MAT *>(&A.mat()), A.ai() + ai, A.aj() + aj,
                        const_cast<VEC *>(&x.vec()), x.ai() + xi, &z.vec(), z.ai() + zi);
    }

    // void fatrop_lu_fact_transposed(const Index m, const Index n, const Index n_max, Index &rank,
    // MAT *At, PermutationMatrix &Pl, PermutationMatrix &Pr, double tol = 1e-5);

    static inline void lu_fact_transposed(const Index m, const Index n, const Index n_max,
                                          Index &rank, const MatRealView &At, PermutationMatrix &Pl,
                                          PermutationMatrix &Pr, double tol = 1e-5)
    {
        fatrop_lu_fact_transposed(m, n, n_max, rank, const_cast<MAT *>(&At.mat()), Pl, Pr, tol);
    }
}

#endif // __fatrop_linear_algebra_blasfeo_operations_hpp__
