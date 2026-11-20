//
// Copyright (c) Lander Vanroye, KU Leuven
//

/**
 * @file lu_factorization.hpp
 * @brief Defines the PermutationMatrix class for LU factorization operations.
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

        int *data() { return permutation_vector_.data(); }

    private:
        const Index size_;                    ///< The size of the permutation matrix.
        std::vector<int> permutation_vector_; ///< The internal representation of the permutation.
    };

    /**
     * @brief LU factorization with complete pivoting of At (n x m, the transpose of A), in place.
     *
     * Only the first n_max rows of At can be pivots. On exit, with r = rank:
     * At[:r, :r] holds L (unit lower) and U, At[:r, r:m] holds U12,
     * At[n_max:n, :r] holds L21 = At21 * inv(U11) and At[n_max:n, r:m] the Schur complement.
     * Pr holds the row swaps (variables), Pl the column swaps (constraints).
     */
    void lu_fact_transposed(const Index m, const Index n, const Index n_max, Index &rank,
                            MatRealView &At, PermutationMatrix &Pl, PermutationMatrix &Pr,
                            Scalar tol);

} // namespace fatrop

#endif // __fatrop_linear_algebra_lu_factorization_hpp__
