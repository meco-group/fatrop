//
// Copyright (c) Lander Vanroye, KU Leuven
//

#include "fatrop/common/exception.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"

namespace fatrop
{
    PermutationMatrix::PermutationMatrix(const Index size) : size_(size), permutation_vector_(size)
    {
        for (Index i = 0; i < size; i++)
        {
            permutation_vector_[i] = i;
        }
    };
    void PermutationMatrix::apply_on_rows(const Index kmax, MAT *mat)
    {
        fatrop_dbg_assert(kmax <= size_);
        ROWPE(kmax, permutation_vector_.data(), mat->n, mat, 0, 0);
    };
    void PermutationMatrix::apply_on_cols(const Index kmax, MAT *mat)
    {
        fatrop_dbg_assert(kmax <= size_);
        COLPE(kmax, permutation_vector_.data(), mat->m, mat, 0, 0);
    };
    void PermutationMatrix::apply(const Index kmax, VEC *vec, const Index ai)
    {
        fatrop_dbg_assert(kmax <= size_);
        fatrop_dbg_assert(ai >= 0);
        VECPE(kmax, permutation_vector_.data(), vec, ai);
    };
    void PermutationMatrix::apply_inverse(const Index kmax, VEC *vec, const Index ai)
    {
        fatrop_dbg_assert(kmax <= size_);
        fatrop_dbg_assert(ai >= 0);
        VECPEI(kmax, permutation_vector_.data(), vec, ai);
    };
    void lu_fact_transposed(const Index m, const Index n, const Index n_max, Index &rank,
                            MatRealView &At, PermutationMatrix &Pl, PermutationMatrix &Pr,
                            Scalar tol)
    {
        fatrop_dbg_assert(n_max <= n && n_max <= At.m() && m <= At.n());
        blasfeo_getrf_fp_wrap(n_max, m, tol, &At.mat(), At.ai(), At.aj(), &At.mat(), At.ai(),
                              At.aj(), Pr.data(), Pl.data(), &rank);
        if (rank == 0)
            return;
        // the rows below n_max were not touched by getrf_fp
        blasfeo_colpe_wrap(rank, Pl.data(), n - n_max, &At.mat(), At.ai() + n_max, At.aj());
        trsm_runn(n - n_max, rank, 1.0, At, 0, 0, At, n_max, 0, At, n_max, 0);
        if (m > rank)
            gemm_nn(n - n_max, m - rank, rank, -1.0, At, n_max, 0, At, 0, rank, 1.0, At, n_max,
                    rank, At, n_max, rank);
    }
}
