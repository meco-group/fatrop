//
// Copyright (c) Lander Vanroye, KU Leuven
//

/**
 * @file lu_factorization.hpp
 * @brief Defines the PermutationMatrix class for LU factorization operations.
 *
 * This file also implements various linear algebra operations required by the FATROP algorithm
 * that are not currently available in the BLASFEO library.
 */

#ifndef __fatrop_linear_algebra_lu_factorization_hpp__
#define __fatrop_linear_algebra_lu_factorization_hpp__

#include "fatrop/common/exception.hpp"
#include "fatrop/context/context.hpp"
#include "fwd.hpp"
#include <vector>

namespace fatrop
{
    /**
     * @class PermutationMatrix
     * @brief Represents a permutation matrix for use in LU factorization.
     *
     * This class provides methods to apply row and column permutations
     * to matrices and vectors, which is essential for LU factorization
     * algorithms to maintain numerical stability and efficiency.
     */
    class PermutationMatrix
    {
    public:
        /**
         * @brief Constructs a PermutationMatrix of the given size.
         * @param size The size of the permutation matrix (number of rows/columns).
         */
        PermutationMatrix(const Index size);

        /**
         * @brief Applies the permutation to the rows of the given matrix.
         * @param kmax The number of permutations to be performed.
         * @param mat Pointer to the matrix to be permuted.
         */
        void apply_on_rows(const Index kmax, MAT *mat);

        /**
         * @brief Applies the permutation to the columns of the given matrix.
         * @param kmax The number of permutations to be performed.
         * @param mat Pointer to the matrix to be permuted.
         */
        void apply_on_cols(const Index kmax, MAT *mat);

        /**
         * @brief Applies the permutation to the given vector.
         * @param kmax The number of permutations to be performed.
         * @param vec Pointer to the vector to be permuted.
         * @param ai The starting index for applying the permutation.
         */
        void apply(const Index kmax, VEC *vec, const Index ai);

        /**
         * @brief Applies the inverse of the permutation to the given vector.
         * @param kmax The number of permutations to be performed.
         * @param vec Pointer to the vector to be permuted.
         * @param ai The starting index for applying the inverse permutation.
         */
        void apply_inverse(const Index kmax, VEC *vec, const Index ai);

        int &operator[](const Index i)
        {
            fatrop_dbg_assert(i >= 0 && i < size_);
            return permutation_vector_[i];
        }

    private:
        const Index size_;                    ///< The size of the permutation matrix.
        std::vector<int> permutation_vector_; ///< The internal representation of the permutation.
    };

    /**
     * @brief Computes LU factorization of a matrix stored in transposed form.
     *
     * Performs LU factorization on a matrix \( A \), where the input \( At \) is \( A^T \)
     * (transposed). Results (\( L \) and \( U \)) are also stored in transposed form. Indices and
     * dimensions refer to the original matrix \( A \).
     *
     * @param m        [in]  Rows of the original matrix \( A \).
     * @param n        [in]  Columns of the original matrix \( A \).
     * @param n_max    [in]  Maximum column storage capacity.
     * @param rank     [out] Effective rank of \( A \).
     * @param At       [in]  Transposed input matrix \( A^T \).
     * @param Pl_p     [out] Left permutation matrix \( P_L \).
     * @param Pr_p     [out] Right permutation matrix \( P_R \).
     * @param tol      [in]  Tolerance for rank determination (default \( 1 \times 10^{-5} \)).
     *
     * @note Ensure \( At \), \( Pl_p \), and \( Pr_p \) are allocated before use.
     */
    void LU_FACT_transposed(const Index m, const Index n, const Index n_max, Index &rank, MAT *At,
                            PermutationMatrix &Pl, PermutationMatrix &Pr, double tol = 1e-5);

    /**
     * @brief Performs an addition operation with a transposed matrix.
     *
     * @param m Number of rows in the result matrix.
     * @param n Number of columns in the result matrix.
     * @param alpha Scalar multiplier for the addition.
     * @param sA Source matrix A.
     * @param offs_ai Row offset for matrix A.
     * @param offs_aj Column offset for matrix A.
     * @param sB Destination matrix B.
     * @param offs_bi Row offset for matrix B.
     * @param offs_bj Column offset for matrix B.
     */
    void fatrop_gead_transposed(Index m, Index n, Scalar alpha, struct blasfeo_dmat *sA,
                                Index offs_ai, Index offs_aj, struct blasfeo_dmat *sB,
                                Index offs_bi, Index offs_bj);

    /**
     * @brief Solves a triangular system of linear equations (upper triangular, not transposed, unit
     * diagonal).
     *
     */
    void fatrop_dtrsv_unu(const Index m, const Index n, blasfeo_dmat *sA, const Index ai,
                          const Index aj, blasfeo_dvec *sx, const Index xi, blasfeo_dvec *sz,
                          const Index zi);

    /**
     * @brief Solves a triangular system of linear equations (upper triangular, transposed, unit
     * diagonal).
     *
     */
    void fatrop_dtrsv_utu(const Index m, blasfeo_dmat *sA, const Index ai, const Index aj,
                          blasfeo_dvec *sx, const Index xi, blasfeo_dvec *sz, const Index zi);

} // namespace fatrop

#endif // __fatrop_linear_algebra_lu_factorization_hpp__
