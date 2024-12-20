//
// Copyright (c) Lander Vanroye, KU Leuven
//

#include "fatrop/linear_algebra/lu_factorization.hpp"
#include "fatrop/common/exception.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include <algorithm>
#include <utility>

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
        ROWPE(kmax, permutation_vector_.data(), mat);
    };
    void PermutationMatrix::apply_on_cols(const Index kmax, MAT *mat)
    {
        fatrop_dbg_assert(kmax <= size_);
        COLPE(kmax, permutation_vector_.data(), mat);
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

    /**
     * \brief Finds the maximum element (by absolute value) in a BLASFEO matrix within a submatrix.
     *
     * This function searches for the maximum absolute value element in a submatrix of size (m, n)
     * of the BLASFEO matrix, starting from the position (ai, aj). It returns the indices of the
     * maximum element as a pair (row, column).
     */
    std::pair<Index, Index> max_el(const Index m, const Index n, const MAT *matr, const Index ai,
                                   const Index aj)
    {
        // Initialize result indices to the starting position
        std::pair<Index, Index> max_indices{ai, aj};

        // Set the initial maximum value to the absolute value at the starting position
        double max_value = abs(MATEL(matr, ai, aj));

        // Iterate over the submatrix to find the maximum element
        for (Index col = aj; col < n; ++col)
        {
            for (Index row = ai; row < m; ++row)
            {
                double current_value = abs(MATEL(matr, row, col));
                if (current_value > max_value)
                {
                    max_value = current_value;
                    max_indices.first = row;
                    max_indices.second = col;
                }
            }
        }
        return max_indices;
    }
    void LU_FACT_transposed(const Index m, const Index n, const Index n_max, Index &rank, MAT *At,
                            PermutationMatrix &Pl, PermutationMatrix &Pr, double tol)
    {
        At->use_dA = 0;
        Index minmn = std::min(m, n_max);
        Index j = 0;
        for (Index i = 0; i < minmn; i++)
        {
            std::pair<Index, Index> max_curr = max_el(n_max, m, At, i, i);
            if (abs(MATEL(At, max_curr.first, max_curr.second)) < tol)
            {
                break;
            }
            // switch rows
            COLSW(n, At, 0, i, At, 0, max_curr.second);
            // save in permutation vector
            Pl[i] = max_curr.second;
            // switch cols
            ROWSW(m, At, i, 0, At, max_curr.first, 0);
            // save in permutation vector
            Pr[i] = max_curr.first;
            for (Index j = i + 1; j < m; j++)
            {
                double Lji = MATEL(At, i, j) / MATEL(At, i, i);
                MATEL(At, i, j) = Lji;
                GEAD(n - (i + 1), 1, -Lji, At, i + 1, i, At, i + 1, j);
            }
            j = i + 1;
        }
        rank = j;
    }
    // B <= B + alpha*A^T (B is mxn)
    void fatrop_gead_transposed(Index m, Index n, Scalar alpha, struct blasfeo_dmat *sA,
                                Index offs_ai, Index offs_aj, struct blasfeo_dmat *sB,
                                Index offs_bi, Index offs_bj)
    {
        for (Index bj = 0; bj < n; bj++)
        {
            for (Index bi = 0; bi < m; bi++)
            {
                MATEL(sB, offs_bi + bi, offs_bj + bj) +=
                    alpha * MATEL(sA, offs_ai + bj, offs_aj + bi);
            }
        }
    }
    void fatrop_dtrsv_unu(const Index m, const Index n, blasfeo_dmat *sA, const Index ai,
                          const Index aj, blasfeo_dvec *sx, const Index xi, blasfeo_dvec *sz,
                          const Index zi)
    {
        for (Index i = m; i < n; i++)
        {
            VECEL(sz, zi + i) = VECEL(sx, xi + i);
        }
        for (Index i = m - 1; i >= 0; i--)
        {
            Scalar res = VECEL(sx, xi + i);
            for (Index j = i + 1; j < n; j++)
            {
                res -= MATEL(sA, ai + i, aj + j) * VECEL(sz, zi + j);
            }
            VECEL(sz, zi + i) = res;
        }
    }

    void fatrop_dtrsv_utu(const Index m, blasfeo_dmat *sA, const Index ai, const Index aj,
                          blasfeo_dvec *sx, const Index xi, blasfeo_dvec *sz, const Index zi)
    {
        for (Index i = 0; i < m; i++)
        {
            Scalar res = VECEL(sx, xi + i);
            for (Index j = 0; j < i; j++)
            {
                res -= MATEL(sA, ai + j, aj + i) * VECEL(sz, zi + j);
            }
            VECEL(sz, zi + i) = res;
        }
    }
}
