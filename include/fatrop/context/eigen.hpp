//
// Copyright (c) 2026 fatrop contributors
//
#ifndef __fatrop_context_eigen_hpp__
#define __fatrop_context_eigen_hpp__

/**
 * @file eigen.hpp
 * @brief Eigen-based linear-algebra backend.
 *
 * Provides plain column-major storage types and implementations of the kernel
 * set fatrop uses, with the same names, signatures and semantics as the
 * BLASFEO kernels mapped in blasfeo.hpp (see that file for the macro list).
 * The storage structs deliberately mirror the blasfeo_dmat/dvec field names
 * (m, n, pA/pa, mem, memsize, use_dA) so code written against those fields
 * compiles unchanged.
 */

#include "fatrop/common/exception.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>

// elementary types
namespace fatrop
{
    typedef double Scalar;
    typedef int Index;
} // namespace fatrop

// column-major storage types (global scope, mirroring blasfeo struct layout conventions)
struct fatrop_eigen_mat
{
    int m;        ///< number of rows (also the leading dimension of pA)
    int n;        ///< number of columns
    double *pA;   ///< column-major data, ld == m
    int use_dA;   ///< unused; present for source compatibility with blasfeo_dmat
    void *mem;    ///< owned allocation (freed by FREE_MAT)
    int memsize;  ///< allocation size in bytes
};

struct fatrop_eigen_vec
{
    int m;        ///< number of elements
    double *pa;   ///< data
    void *mem;    ///< owned allocation (freed by FREE_VEC)
    int memsize;  ///< allocation size in bytes
};

namespace fatrop_eigen
{
    typedef Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, Eigen::OuterStride<>> MapMat;
    typedef Eigen::Map<const Eigen::MatrixXd, Eigen::Unaligned, Eigen::OuterStride<>> CMapMat;
    typedef Eigen::Map<Eigen::VectorXd> MapVec;
    typedef Eigen::Map<const Eigen::VectorXd> CMapVec;

    inline MapMat blk(fatrop_eigen_mat *A, int ai, int aj, int m, int n)
    {
        return MapMat(A->pA + static_cast<std::ptrdiff_t>(aj) * A->m + ai, m, n,
                      Eigen::OuterStride<>(A->m));
    }
    inline CMapMat blk(const fatrop_eigen_mat *A, int ai, int aj, int m, int n)
    {
        return CMapMat(A->pA + static_cast<std::ptrdiff_t>(aj) * A->m + ai, m, n,
                       Eigen::OuterStride<>(A->m));
    }
    inline MapVec seg(fatrop_eigen_vec *x, int xi, int m) { return MapVec(x->pa + xi, m); }
    inline CMapVec seg(const fatrop_eigen_vec *x, int xi, int m)
    {
        return CMapVec(x->pa + xi, m);
    }

    // grow-only thread-local scratch buffer for kernels that need a temporary
    inline double *scratch(std::size_t nelem)
    {
        static thread_local std::vector<double> buf;
        if (buf.size() < nelem)
            buf.resize(nelem);
        return buf.data();
    }
} // namespace fatrop_eigen

// ---------------------------------------------------------------------------
// allocation
// ---------------------------------------------------------------------------

inline int fatrop_eigen_memsize_mat(int m, int n)
{
    size_t bytes = static_cast<size_t>(m) * static_cast<size_t>(n) * sizeof(double);
    bytes = (bytes + 63) / 64 * 64;
    if (bytes == 0)
        bytes = 64;
    return static_cast<int>(bytes);
}

inline int fatrop_eigen_memsize_vec(int m)
{
    size_t bytes = static_cast<size_t>(m) * sizeof(double);
    bytes = (bytes + 63) / 64 * 64;
    if (bytes == 0)
        bytes = 64;
    return static_cast<int>(bytes);
}

inline void fatrop_eigen_allocate_mat(int m, int n, fatrop_eigen_mat *sA)
{
    sA->m = m;
    sA->n = n;
    sA->use_dA = 0;
    sA->memsize = fatrop_eigen_memsize_mat(m, n);
    sA->mem = std::aligned_alloc(64, static_cast<size_t>(sA->memsize));
    sA->pA = static_cast<double *>(sA->mem);
}

inline void fatrop_eigen_free_mat(fatrop_eigen_mat *sA) { std::free(sA->mem); }

inline void fatrop_eigen_allocate_vec(int m, fatrop_eigen_vec *sx)
{
    sx->m = m;
    sx->memsize = fatrop_eigen_memsize_vec(m);
    sx->mem = std::aligned_alloc(64, static_cast<size_t>(sx->memsize));
    sx->pa = static_cast<double *>(sx->mem);
}

inline void fatrop_eigen_free_vec(fatrop_eigen_vec *sx) { std::free(sx->mem); }

// ---------------------------------------------------------------------------
// element access
// ---------------------------------------------------------------------------

inline double &fatrop_eigen_matel(fatrop_eigen_mat *A, int ai, int aj)
{
    return A->pA[static_cast<std::ptrdiff_t>(aj) * A->m + ai];
}
inline double fatrop_eigen_matel(const fatrop_eigen_mat *A, int ai, int aj)
{
    return A->pA[static_cast<std::ptrdiff_t>(aj) * A->m + ai];
}
inline double &fatrop_eigen_vecel(fatrop_eigen_vec *x, int xi) { return x->pa[xi]; }
inline double fatrop_eigen_vecel(const fatrop_eigen_vec *x, int xi) { return x->pa[xi]; }

// ---------------------------------------------------------------------------
// vector kernels
// ---------------------------------------------------------------------------

// z = alpha * x + beta * y
inline void fatrop_eigen_axpby(int m, double alpha, fatrop_eigen_vec *x, int xi, double beta,
                               fatrop_eigen_vec *y, int yi, fatrop_eigen_vec *z, int zi)
{
    fatrop_eigen::seg(z, zi, m) =
        alpha * fatrop_eigen::seg(x, xi, m) + beta * fatrop_eigen::seg(y, yi, m);
}

// z = alpha * x + y
inline void fatrop_eigen_axpy(int m, double alpha, fatrop_eigen_vec *x, int xi,
                              fatrop_eigen_vec *y, int yi, fatrop_eigen_vec *z, int zi)
{
    fatrop_eigen::seg(z, zi, m) = alpha * fatrop_eigen::seg(x, xi, m) + fatrop_eigen::seg(y, yi, m);
}

inline void fatrop_eigen_vecse(int m, double alpha, fatrop_eigen_vec *x, int xi)
{
    fatrop_eigen::seg(x, xi, m).setConstant(alpha);
}

inline void fatrop_eigen_vecsc(int m, double alpha, fatrop_eigen_vec *x, int xi)
{
    fatrop_eigen::seg(x, xi, m) *= alpha;
}

// y = alpha * x
inline void fatrop_eigen_veccpsc(int m, double alpha, fatrop_eigen_vec *x, int xi,
                                 fatrop_eigen_vec *y, int yi)
{
    fatrop_eigen::seg(y, yi, m) = alpha * fatrop_eigen::seg(x, xi, m);
}

inline void fatrop_eigen_veccp(int m, fatrop_eigen_vec *x, int xi, fatrop_eigen_vec *y, int yi)
{
    fatrop_eigen::seg(y, yi, m) = fatrop_eigen::seg(x, xi, m);
}

// z = x .* y
inline void fatrop_eigen_vecmul(int m, fatrop_eigen_vec *x, int xi, fatrop_eigen_vec *y, int yi,
                                fatrop_eigen_vec *z, int zi)
{
    fatrop_eigen::seg(z, zi, m) =
        fatrop_eigen::seg(x, xi, m).cwiseProduct(fatrop_eigen::seg(y, yi, m));
}

// z += x .* y
inline void fatrop_eigen_vecmulacc(int m, fatrop_eigen_vec *x, int xi, fatrop_eigen_vec *y,
                                   int yi, fatrop_eigen_vec *z, int zi)
{
    fatrop_eigen::seg(z, zi, m) +=
        fatrop_eigen::seg(x, xi, m).cwiseProduct(fatrop_eigen::seg(y, yi, m));
}

inline double fatrop_eigen_dot(int m, fatrop_eigen_vec *x, int xi, fatrop_eigen_vec *y, int yi)
{
    return fatrop_eigen::seg(x, xi, m).dot(fatrop_eigen::seg(y, yi, m));
}

// forward sequential element swaps (LAPACK-style)
inline void fatrop_eigen_vecpe(int kmax, int *ipiv, fatrop_eigen_vec *sx, int xi)
{
    double *x = sx->pa + xi;
    for (int i = 0; i < kmax; ++i)
        if (ipiv[i] != i)
            std::swap(x[i], x[ipiv[i]]);
}

inline void fatrop_eigen_vecpei(int kmax, int *ipiv, fatrop_eigen_vec *sx, int xi)
{
    double *x = sx->pa + xi;
    for (int i = kmax - 1; i >= 0; --i)
        if (ipiv[i] != i)
            std::swap(x[i], x[ipiv[i]]);
}

inline void fatrop_eigen_pack_vec(int m, double *x, int incx, fatrop_eigen_vec *sx, int xi)
{
    for (int i = 0; i < m; ++i)
        sx->pa[xi + i] = x[i * incx];
}

inline void fatrop_eigen_unpack_vec(int m, fatrop_eigen_vec *sx, int xi, double *x, int incx)
{
    for (int i = 0; i < m; ++i)
        x[i * incx] = sx->pa[xi + i];
}

// ---------------------------------------------------------------------------
// matrix data-movement kernels
// ---------------------------------------------------------------------------

inline void fatrop_eigen_gese(int m, int n, double alpha, fatrop_eigen_mat *sA, int ai, int aj)
{
    fatrop_eigen::blk(sA, ai, aj, m, n).setConstant(alpha);
}

inline void fatrop_eigen_gecp(int m, int n, fatrop_eigen_mat *sA, int ai, int aj,
                              fatrop_eigen_mat *sB, int bi, int bj)
{
    fatrop_eigen::blk(sB, bi, bj, m, n) = fatrop_eigen::blk(sA, ai, aj, m, n);
}

// B += alpha * A
inline void fatrop_eigen_gead(int m, int n, double alpha, fatrop_eigen_mat *sA, int ai, int aj,
                              fatrop_eigen_mat *sB, int bi, int bj)
{
    fatrop_eigen::blk(sB, bi, bj, m, n) += alpha * fatrop_eigen::blk(sA, ai, aj, m, n);
}

inline void fatrop_eigen_gesc(int m, int n, double alpha, fatrop_eigen_mat *sA, int ai, int aj)
{
    fatrop_eigen::blk(sA, ai, aj, m, n) *= alpha;
}

// B = A^T (no aliasing support, as in BLASFEO)
inline void fatrop_eigen_getr(int m, int n, fatrop_eigen_mat *sA, int ai, int aj,
                              fatrop_eigen_mat *sB, int bi, int bj)
{
    fatrop_eigen::blk(sB, bi, bj, n, m) = fatrop_eigen::blk(sA, ai, aj, m, n).transpose();
}

// B(j,i) = A(i,j) for i >= j: transpose of the lower triangle of A written into
// the upper triangle of B, other entries of B untouched (BLASFEO dtrtr_l).
inline void fatrop_eigen_trtr_l(int m, fatrop_eigen_mat *sA, int ai, int aj,
                                fatrop_eigen_mat *sB, int bi, int bj)
{
    for (int j = 0; j < m; ++j)
        for (int i = j; i < m; ++i)
            fatrop_eigen_matel(sB, bi + j, bj + i) = fatrop_eigen_matel(sA, ai + i, aj + j);
}

// A(i,i) += alpha (BLASFEO ddiare adds)
inline void fatrop_eigen_diare(int kmax, double alpha, fatrop_eigen_mat *sA, int ai, int aj)
{
    for (int i = 0; i < kmax; ++i)
        fatrop_eigen_matel(sA, ai + i, aj + i) += alpha;
}

// A(i,i) += alpha * x[i]
inline void fatrop_eigen_diaad(int kmax, double alpha, fatrop_eigen_vec *sx, int xi,
                               fatrop_eigen_mat *sA, int ai, int aj)
{
    for (int i = 0; i < kmax; ++i)
        fatrop_eigen_matel(sA, ai + i, aj + i) += alpha * sx->pa[xi + i];
}

inline void fatrop_eigen_rowex(int kmax, double alpha, fatrop_eigen_mat *sA, int ai, int aj,
                               fatrop_eigen_vec *sx, int xi)
{
    for (int i = 0; i < kmax; ++i)
        sx->pa[xi + i] = alpha * fatrop_eigen_matel(sA, ai, aj + i);
}

inline void fatrop_eigen_rowin(int kmax, double alpha, fatrop_eigen_vec *sx, int xi,
                               fatrop_eigen_mat *sA, int ai, int aj)
{
    for (int i = 0; i < kmax; ++i)
        fatrop_eigen_matel(sA, ai, aj + i) = alpha * sx->pa[xi + i];
}

inline void fatrop_eigen_colin(int kmax, fatrop_eigen_vec *sx, int xi, fatrop_eigen_mat *sA,
                               int ai, int aj)
{
    fatrop_eigen::blk(sA, ai, aj, kmax, 1).col(0) = fatrop_eigen::seg(sx, xi, kmax);
}

inline void fatrop_eigen_colsc(int kmax, double alpha, fatrop_eigen_mat *sA, int ai, int aj)
{
    fatrop_eigen::blk(sA, ai, aj, kmax, 1) *= alpha;
}

// swap kmax elements of row ai (from col aj) of A with row ci (from col cj) of C
inline void fatrop_eigen_rowsw(int kmax, fatrop_eigen_mat *sA, int ai, int aj,
                               fatrop_eigen_mat *sC, int ci, int cj)
{
    for (int i = 0; i < kmax; ++i)
        std::swap(fatrop_eigen_matel(sA, ai, aj + i), fatrop_eigen_matel(sC, ci, cj + i));
}

// swap kmax elements of col aj (from row ai) of A with col cj (from row ci) of C
inline void fatrop_eigen_colsw(int kmax, fatrop_eigen_mat *sA, int ai, int aj,
                               fatrop_eigen_mat *sC, int ci, int cj)
{
    for (int i = 0; i < kmax; ++i)
        std::swap(fatrop_eigen_matel(sA, ai + i, aj), fatrop_eigen_matel(sC, ci + i, cj));
}

// LAPACK-style forward sequential full-row swaps
inline void fatrop_eigen_rowpe(int kmax, int *ipiv, fatrop_eigen_mat *sA)
{
    for (int i = 0; i < kmax; ++i)
        if (ipiv[i] != i)
            fatrop_eigen_rowsw(sA->n, sA, i, 0, sA, ipiv[i], 0);
}

inline void fatrop_eigen_rowpei(int kmax, int *ipiv, fatrop_eigen_mat *sA)
{
    for (int i = kmax - 1; i >= 0; --i)
        if (ipiv[i] != i)
            fatrop_eigen_rowsw(sA->n, sA, i, 0, sA, ipiv[i], 0);
}

inline void fatrop_eigen_colpe(int kmax, int *ipiv, fatrop_eigen_mat *sA)
{
    for (int i = 0; i < kmax; ++i)
        if (ipiv[i] != i)
            fatrop_eigen_colsw(sA->m, sA, 0, i, sA, 0, ipiv[i]);
}

inline void fatrop_eigen_colpei(int kmax, int *ipiv, fatrop_eigen_mat *sA)
{
    for (int i = kmax - 1; i >= 0; --i)
        if (ipiv[i] != i)
            fatrop_eigen_colsw(sA->m, sA, 0, i, sA, 0, ipiv[i]);
}

inline void fatrop_eigen_pack_mat(int m, int n, double *A, int lda, fatrop_eigen_mat *sA, int ai,
                                  int aj)
{
    fatrop_eigen::blk(sA, ai, aj, m, n) =
        Eigen::Map<const Eigen::MatrixXd, Eigen::Unaligned, Eigen::OuterStride<>>(
            A, m, n, Eigen::OuterStride<>(lda));
}

inline void fatrop_eigen_unpack_mat(int m, int n, fatrop_eigen_mat *sA, int ai, int aj, double *A,
                                    int lda)
{
    Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, Eigen::OuterStride<>>(
        A, m, n, Eigen::OuterStride<>(lda)) = fatrop_eigen::blk(sA, ai, aj, m, n);
}

// ---------------------------------------------------------------------------
// BLAS-level kernels
// ---------------------------------------------------------------------------

// D = beta * C + alpha * A * B^T   (A m x k, B n x k)
inline void fatrop_eigen_gemm_nt(int m, int n, int k, double alpha, fatrop_eigen_mat *sA, int ai,
                                 int aj, fatrop_eigen_mat *sB, int bi, int bj, double beta,
                                 fatrop_eigen_mat *sC, int ci, int cj, fatrop_eigen_mat *sD,
                                 int di, int dj)
{
    if (m <= 0 || n <= 0)
        return;
    auto A = fatrop_eigen::blk(sA, ai, aj, m, k);
    auto B = fatrop_eigen::blk(sB, bi, bj, n, k);
    auto D = fatrop_eigen::blk(sD, di, dj, m, n);
    if (beta == 0.0)
    {
        // C must not be read when beta == 0 (may be uninitialized)
        D.noalias() = alpha * (A * B.transpose());
    }
    else
    {
        auto C = fatrop_eigen::blk(sC, ci, cj, m, n);
        fatrop_eigen::MapMat T(fatrop_eigen::scratch(static_cast<std::size_t>(m) * n), m, n,
                               Eigen::OuterStride<>(m));
        T.noalias() = alpha * (A * B.transpose());
        D = beta * C + T; // elementwise; safe when D aliases C
    }
}

// lower-trapezoid syrk: D(i,j) = beta*C(i,j) + alpha*(A*B^T)(i,j) for 0<=j<n, j<=i<m
inline void fatrop_eigen_syrk_ln_mn(int m, int n, int k, double alpha, fatrop_eigen_mat *sA,
                                    int ai, int aj, fatrop_eigen_mat *sB, int bi, int bj,
                                    double beta, fatrop_eigen_mat *sC, int ci, int cj,
                                    fatrop_eigen_mat *sD, int di, int dj)
{
    if (m <= 0 || n <= 0)
        return;
    auto A = fatrop_eigen::blk(sA, ai, aj, m, k);
    auto B = fatrop_eigen::blk(sB, bi, bj, n, k);
    auto D = fatrop_eigen::blk(sD, di, dj, m, n);
    fatrop_eigen::MapMat T(fatrop_eigen::scratch(static_cast<std::size_t>(m) * n), m, n,
                           Eigen::OuterStride<>(m));
    T.noalias() = alpha * (A * B.transpose());
    if (beta == 0.0)
    {
        for (int j = 0; j < n && j < m; ++j)
            D.col(j).segment(j, m - j) = T.col(j).segment(j, m - j);
    }
    else
    {
        auto C = fatrop_eigen::blk(sC, ci, cj, m, n);
        for (int j = 0; j < n && j < m; ++j)
            D.col(j).segment(j, m - j) =
                beta * C.col(j).segment(j, m - j) + T.col(j).segment(j, m - j);
    }
}

inline void fatrop_eigen_syrk_ln(int m, int k, double alpha, fatrop_eigen_mat *sA, int ai, int aj,
                                 fatrop_eigen_mat *sB, int bi, int bj, double beta,
                                 fatrop_eigen_mat *sC, int ci, int cj, fatrop_eigen_mat *sD,
                                 int di, int dj)
{
    fatrop_eigen_syrk_ln_mn(m, m, k, alpha, sA, ai, aj, sB, bi, bj, beta, sC, ci, cj, sD, di, dj);
}

// Cholesky of a lower m x n trapezoid (m >= n): factor the top n x n block,
// triangular-solve the rows below. Reproduces the BLASFEO convention of writing
// a zero column when a pivot is not strictly positive (callers detect this via
// a small-diagonal check instead of an error code).
inline void fatrop_eigen_potrf_l_mn(int m, int n, fatrop_eigen_mat *sC, int ci, int cj,
                                    fatrop_eigen_mat *sD, int di, int dj)
{
    if (m <= 0 || n <= 0)
        return;
    // copy the lower trapezoid of C into D if they differ
    if (sD->pA != sC->pA || di != ci || dj != cj)
    {
        for (int j = 0; j < n && j < m; ++j)
            fatrop_eigen::blk(sD, di + j, dj + j, m - j, 1) =
                fatrop_eigen::blk(sC, ci + j, cj + j, m - j, 1);
    }
    auto D = fatrop_eigen::blk(sD, di, dj, m, n);
    for (int j = 0; j < n; ++j)
    {
        double d = D(j, j) - D.row(j).head(j).squaredNorm();
        const int mm = m - j - 1;
        if (d > 0.0)
        {
            const double l = std::sqrt(d);
            D(j, j) = l;
            if (mm > 0)
            {
                auto col = D.col(j).segment(j + 1, mm);
                if (j > 0)
                    col.noalias() -= D.block(j + 1, 0, mm, j) * D.row(j).head(j).transpose();
                col /= l;
            }
        }
        else
        {
            D.col(j).segment(j, mm + 1).setZero();
        }
    }
}

// D = alpha * B * A^{-T}, A n x n lower triangular (non-unit diagonal)
inline void fatrop_eigen_trsm_rltn(int m, int n, double alpha, fatrop_eigen_mat *sA, int ai,
                                   int aj, fatrop_eigen_mat *sB, int bi, int bj,
                                   fatrop_eigen_mat *sD, int di, int dj)
{
    if (m <= 0 || n <= 0)
        return;
    auto A = fatrop_eigen::blk(sA, ai, aj, n, n);
    auto B = fatrop_eigen::blk(sB, bi, bj, m, n);
    auto D = fatrop_eigen::blk(sD, di, dj, m, n);
    D = alpha * B; // safe when D aliases B
    A.transpose().triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheRight>(D);
}

// D = alpha * B * A^{-1}, A n x n lower triangular (non-unit diagonal)
inline void fatrop_eigen_trsm_rlnn(int m, int n, double alpha, fatrop_eigen_mat *sA, int ai,
                                   int aj, fatrop_eigen_mat *sB, int bi, int bj,
                                   fatrop_eigen_mat *sD, int di, int dj)
{
    if (m <= 0 || n <= 0)
        return;
    auto A = fatrop_eigen::blk(sA, ai, aj, n, n);
    auto B = fatrop_eigen::blk(sB, bi, bj, m, n);
    auto D = fatrop_eigen::blk(sD, di, dj, m, n);
    D = alpha * B;
    A.triangularView<Eigen::Lower>().solveInPlace<Eigen::OnTheRight>(D);
}

// z = A^{-1} x, A lower triangular
inline void fatrop_eigen_trsv_lnn(int m, fatrop_eigen_mat *sA, int ai, int aj,
                                  fatrop_eigen_vec *sx, int xi, fatrop_eigen_vec *sz, int zi)
{
    if (m <= 0)
        return;
    auto A = fatrop_eigen::blk(sA, ai, aj, m, m);
    auto z = fatrop_eigen::seg(sz, zi, m);
    z = fatrop_eigen::seg(sx, xi, m);
    A.triangularView<Eigen::Lower>().solveInPlace(z);
}

// z = A^{-T} x, A lower triangular
inline void fatrop_eigen_trsv_ltn(int m, fatrop_eigen_mat *sA, int ai, int aj,
                                  fatrop_eigen_vec *sx, int xi, fatrop_eigen_vec *sz, int zi)
{
    if (m <= 0)
        return;
    auto A = fatrop_eigen::blk(sA, ai, aj, m, m);
    auto z = fatrop_eigen::seg(sz, zi, m);
    z = fatrop_eigen::seg(sx, xi, m);
    A.transpose().triangularView<Eigen::Upper>().solveInPlace(z);
}

// z = A^{-T} x, A upper triangular
inline void fatrop_eigen_trsv_utn(int m, fatrop_eigen_mat *sA, int ai, int aj,
                                  fatrop_eigen_vec *sx, int xi, fatrop_eigen_vec *sz, int zi)
{
    if (m <= 0)
        return;
    auto A = fatrop_eigen::blk(sA, ai, aj, m, m);
    auto z = fatrop_eigen::seg(sz, zi, m);
    z = fatrop_eigen::seg(sx, xi, m);
    A.transpose().triangularView<Eigen::Lower>().solveInPlace(z);
}

// z = beta * y + alpha * A^T * x   (A m x n, x m, y/z n)
inline void fatrop_eigen_gemv_t(int m, int n, double alpha, fatrop_eigen_mat *sA, int ai, int aj,
                                fatrop_eigen_vec *sx, int xi, double beta, fatrop_eigen_vec *sy,
                                int yi, fatrop_eigen_vec *sz, int zi)
{
    if (n <= 0)
        return;
    auto A = fatrop_eigen::blk(sA, ai, aj, m, n);
    auto x = fatrop_eigen::seg(sx, xi, m);
    auto z = fatrop_eigen::seg(sz, zi, n);
    if (beta == 0.0)
        z = alpha * (A.transpose() * x);
    else
        z = beta * fatrop_eigen::seg(sy, yi, n) + alpha * (A.transpose() * x);
}

// z = beta * y + alpha * A * x   (A m x n, x n, y/z m)
inline void fatrop_eigen_gemv_n(int m, int n, double alpha, fatrop_eigen_mat *sA, int ai, int aj,
                                fatrop_eigen_vec *sx, int xi, double beta, fatrop_eigen_vec *sy,
                                int yi, fatrop_eigen_vec *sz, int zi)
{
    if (m <= 0)
        return;
    auto A = fatrop_eigen::blk(sA, ai, aj, m, n);
    auto x = fatrop_eigen::seg(sx, xi, n);
    auto z = fatrop_eigen::seg(sz, zi, m);
    if (beta == 0.0)
        z = alpha * (A * x);
    else
        z = beta * fatrop_eigen::seg(sy, yi, m) + alpha * (A * x);
}

// D = C + alpha * x * y^T
inline void fatrop_eigen_ger(int m, int n, double alpha, fatrop_eigen_vec *sx, int xi,
                             fatrop_eigen_vec *sy, int yi, fatrop_eigen_mat *sC, int ci, int cj,
                             fatrop_eigen_mat *sD, int di, int dj)
{
    if (m <= 0 || n <= 0)
        return;
    auto D = fatrop_eigen::blk(sD, di, dj, m, n);
    if (sD->pA != sC->pA || di != ci || dj != cj)
        D = fatrop_eigen::blk(sC, ci, cj, m, n);
    D.noalias() +=
        alpha * fatrop_eigen::seg(sx, xi, m) * fatrop_eigen::seg(sy, yi, n).transpose();
}

// ---------------------------------------------------------------------------
// macro mapping (same names as in blasfeo.hpp)
// ---------------------------------------------------------------------------

#define VEC fatrop_eigen_vec
#define VECEL fatrop_eigen_vecel
#define MAT fatrop_eigen_mat
#define MATEL fatrop_eigen_matel

// Vector-related definitions
#define AXPBY fatrop_eigen_axpby
#define AXPY fatrop_eigen_axpy
#define VECPE fatrop_eigen_vecpe
#define VECPEI fatrop_eigen_vecpei
#define VECSE fatrop_eigen_vecse
#define VECSC fatrop_eigen_vecsc
#define VECCPSC fatrop_eigen_veccpsc
#define VECCP fatrop_eigen_veccp
#define VECMUL fatrop_eigen_vecmul
#define DOT fatrop_eigen_dot

#define ALLOCATE_VEC fatrop_eigen_allocate_vec
#define FREE_VEC fatrop_eigen_free_vec
#define MEMSIZE_VEC fatrop_eigen_memsize_vec

// Matrix-related definitions
#define ROWPE fatrop_eigen_rowpe
#define ROWPEI fatrop_eigen_rowpei
#define COLPE fatrop_eigen_colpe
#define COLPEI fatrop_eigen_colpei
#define ROWSW fatrop_eigen_rowsw
#define COLSW fatrop_eigen_colsw
#define GEAD fatrop_eigen_gead
#define GECP fatrop_eigen_gecp
#define GESC fatrop_eigen_gesc
#define TRSM_RLTN fatrop_eigen_trsm_rltn
#define TRSM_RLNN fatrop_eigen_trsm_rlnn
#define GEMM_NT fatrop_eigen_gemm_nt
#define SYRK_LN_MN fatrop_eigen_syrk_ln_mn
#define SYRK_LN fatrop_eigen_syrk_ln
#define GETR fatrop_eigen_getr
#define TRTR_L fatrop_eigen_trtr_l
#define POTRF_L_MN fatrop_eigen_potrf_l_mn
#define ROWEX fatrop_eigen_rowex
#define ROWIN fatrop_eigen_rowin
#define COLIN fatrop_eigen_colin
#define TRSV_LTN fatrop_eigen_trsv_ltn
#define TRSV_LNN fatrop_eigen_trsv_lnn
#define TRSV_UTN fatrop_eigen_trsv_utn
#define GEMV_T fatrop_eigen_gemv_t
#define GEMV_N fatrop_eigen_gemv_n
#define GESE fatrop_eigen_gese
#define DIARE fatrop_eigen_diare
#define DIAAD fatrop_eigen_diaad
#define COLSC fatrop_eigen_colsc
#define VECMULACC fatrop_eigen_vecmulacc
#define GER fatrop_eigen_ger

#define PACK_MAT fatrop_eigen_pack_mat
#define UNPACK_MAT fatrop_eigen_unpack_mat
#define PACK_VEC fatrop_eigen_pack_vec
#define UNPACK_VEC fatrop_eigen_unpack_vec

#define ALLOCATE_MAT fatrop_eigen_allocate_mat
#define FREE_MAT fatrop_eigen_free_mat
#define MEMSIZE_MAT fatrop_eigen_memsize_mat

#endif // __fatrop_context_eigen_hpp__
