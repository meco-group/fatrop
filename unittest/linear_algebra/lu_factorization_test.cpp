#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "../random_matrix.hpp"
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>

namespace fatrop
{
    namespace test
    {
        using fatrop::MatRealAllocated;
        using fatrop::VecRealAllocated;

        TEST(PermutationMatrixTest, Constructor)
        {
            // Test that the constructor initializes the permutation as identity
            const Index size = 5;
            PermutationMatrix perm(size);
            for (Index i = 0; i < size; ++i)
            {
                EXPECT_EQ(perm[i], i);
            }
        }

        TEST(PermutationMatrixTest, ApplyOnRows)
        {
            // Test applying permutation on matrix rows
            const Index size = 3;
            PermutationMatrix perm(size);

            // Create a permutation: swap rows 0 and 2, and rows 2 and 1
            // 0, 1, 2 --> 2, 0, 1
            perm[0] = 2;
            perm[1] = 1;
            perm[2] = 1;
            const Index perm_result[] = {2, 0, 1};

            MatRealAllocated matrix(size, size);
            for (Index i = 0; i < size; ++i)
                for (Index j = 0; j < size; ++j)
                    matrix(i, j) = i * size + j + 1;

            perm.apply_on_rows(size, &matrix.mat());

            VecRealAllocated row(size);
            for (Index i = 0; i < size; ++i)
            {
                row = matrix.row(i);
                for (Index j = 0; j < size; ++j)
                {
                    Scalar expected = perm_result[i] * size + j + 1;
                    EXPECT_DOUBLE_EQ(expected, row(j));
                }
            }
        }

        TEST(PermutationMatrixTest, ApplyOnCols)
        {
            // Test applying permutation on matrix columns
            const Index size = 3;
            PermutationMatrix perm(size);

            // Create a permutation: move column 0 to 2, 1 to 0, and 2 to 1
            perm[0] = 2;
            perm[1] = 0;
            perm[2] = 1;

            MatRealAllocated matrix(size, size);
            for (Index i = 0; i < size; ++i)
                for (Index j = 0; j < size; ++j)
                    matrix(i, j) = i * size + j + 1;

            perm.apply_on_cols(size, &matrix.mat());

            const Index perm_result[] = {1, 0, 2}; // Inverse permutation
            for (Index i = 0; i < size; ++i)
            {
                for (Index j = 0; j < size; ++j)
                {
                    Scalar expected = i * size + perm_result[j] + 1;
                    EXPECT_DOUBLE_EQ(expected, matrix(i, j));
                }
            }
        }

        TEST(PermutationMatrixTest, ApplyOnVector)
        {
            // Test applying permutation on a vector
            const Index size = 3;
            PermutationMatrix perm(size);

            // Create a permutation: move element 0 to 2, 1 to 0, and 2 to 1
            perm[0] = 2;
            perm[1] = 0;
            perm[2] = 1;

            VecRealAllocated vector(size);
            for (Index i = 0; i < size; ++i)
                vector(i) = i + 1;

            perm.apply(size, &vector.vec(), 0);

            const Scalar expected[] = {2, 1, 3};
            for (Index i = 0; i < size; ++i)
            {
                EXPECT_DOUBLE_EQ(expected[i], vector(i));
            }
        }

        TEST(PermutationMatrixTest, ApplyInverseOnVector)
        {
            // Test applying inverse permutation on a vector
            const Index size = 3;
            PermutationMatrix perm(size);

            // Create a permutation: move element 0 to 2, 1 to 0, and 2 to 1
            perm[0] = 2;
            perm[1] = 0;
            perm[2] = 1;

            VecRealAllocated vector(size);
            for (Index i = 0; i < size; ++i)
                vector(i) = i + 1;

            perm.apply(size, &vector.vec(), 0);
            perm.apply_inverse(size, &vector.vec(), 0);

            for (Index i = 0; i < size; ++i)
            {
                EXPECT_DOUBLE_EQ(i + 1, vector(i));
            }
        }

        TEST(PermutationMatrixTest, LargePermutation)
        {
            // Test a larger permutation to ensure it works for non-trivial cases
            const Index size = 5;
            PermutationMatrix perm(size);

            // Create a more complex permutation
            perm[0] = 3;
            perm[1] = 2;
            perm[2] = 4;
            perm[3] = 0;
            perm[4] = 1;

            VecRealAllocated vector(size);
            for (Index i = 0; i < size; ++i)
                vector(i) = i + 1;

            Index perm_result[] = {0, 1, 2, 3, 4};
            for (Index i = 0; i < size; ++i)
            {
                Index tmp = perm_result[i];
                perm_result[i] = perm_result[perm[i]];
                perm_result[perm[i]] = tmp;
            }
            Scalar expected[5];
            for (Index i = 0; i < size; i++)
            {
                expected[i] = vector(perm_result[i]);
            }

            perm.apply(size, &vector.vec(), 0);

            for (Index i = 0; i < size; ++i)
            {
                EXPECT_DOUBLE_EQ(expected[i], vector(i));
            }

            // Apply inverse to get back the original vector
            perm.apply_inverse(size, &vector.vec(), 0);

            for (Index i = 0; i < size; ++i)
            {
                EXPECT_DOUBLE_EQ(i + 1, vector(i));
            }
        }

        // At is n x m, only its first n_max rows may be pivots. Checks Pr * At * Pl = L * U + S,
        // with L (n x r) and U (r x m) read from the factored At and S its trailing block.
        void check_lu_fact_transposed(const Index m, const Index n, const Index n_max,
                                      const MatRealAllocated &At0, const Index rank_expected)
        {
            const Scalar tol = 1e-8;
            MatRealAllocated At(n, m);
            MatRealAllocated Ap(n, m);
            for (Index i = 0; i < n; ++i)
                for (Index j = 0; j < m; ++j)
                    At(i, j) = Ap(i, j) = At0(i, j);
            PermutationMatrix Pl(m);
            PermutationMatrix Pr(n_max);
            Index rank;
            lu_fact_transposed(m, n, n_max, rank, At, Pl, Pr, tol);
            EXPECT_EQ(rank, rank_expected);

            Pr.apply_on_rows(rank, &Ap.mat());
            Pl.apply_on_cols(rank, &Ap.mat());
            for (Index i = 0; i < n; ++i)
            {
                for (Index j = 0; j < m; ++j)
                {
                    Scalar lu = (i >= rank && j >= rank) ? At(i, j) : 0.0;
                    for (Index k = 0; k < std::min({i + 1, j + 1, rank}); ++k)
                    {
                        const Scalar l = (k == i) ? 1.0 : At(i, k);
                        lu += l * At(k, j);
                    }
                    EXPECT_NEAR(Ap(i, j), lu, 1e-10) << "i=" << i << " j=" << j;
                }
            }
            // eligible rows: what is left is below tol
            for (Index i = rank; i < n_max; ++i)
                for (Index j = rank; j < m; ++j)
                    EXPECT_LE(std::abs(At(i, j)), tol);
        }

        MatRealAllocated random_low_rank(Index n, Index m, Index r)
        {
            MatRealAllocated X = random_matrix(n, r, -1.0, 1.0);
            MatRealAllocated Y = random_matrix(r, m, -1.0, 1.0);
            MatRealAllocated A(n, m);
            for (Index i = 0; i < n; ++i)
                for (Index j = 0; j < m; ++j)
                {
                    Scalar sum = 0.0;
                    for (Index k = 0; k < r; ++k)
                        sum += X(i, k) * Y(k, j);
                    A(i, j) = sum;
                }
            return A;
        }

        TEST(LUFactorizationTest, transposedSquare)
        {
            check_lu_fact_transposed(5, 5, 5, random_matrix(5, 5, -1.0, 1.0), 5);
        }

        TEST(LUFactorizationTest, transposedExtraRows)
        {
            // n_max < n: the bottom rows are eliminated but never pivot
            check_lu_fact_transposed(3, 9, 6, random_matrix(9, 3, -1.0, 1.0), 3);
            check_lu_fact_transposed(5, 9, 6, random_matrix(9, 5, -1.0, 1.0), 5);
        }

        TEST(LUFactorizationTest, transposedPivotOnlyEligibleRows)
        {
            // eligible rows tiny, extra rows large: pivots must still come from the top
            MatRealAllocated At = random_matrix(7, 3, -1.0, 1.0);
            for (Index i = 0; i < 4; ++i)
                for (Index j = 0; j < 3; ++j)
                    At(i, j) *= 1e-3;
            check_lu_fact_transposed(3, 7, 4, At, 3);
        }

        TEST(LUFactorizationTest, transposedRankDeficient)
        {
            // rank 2 in the eligible rows, full rank in the extra rows
            MatRealAllocated At = random_low_rank(9, 5, 2);
            for (Index i = 6; i < 9; ++i)
                for (Index j = 0; j < 5; ++j)
                    At(i, j) = random(-1.0, 1.0);
            check_lu_fact_transposed(5, 9, 6, At, 2);
        }

        TEST(LUFactorizationTest, transposedZero)
        {
            MatRealAllocated At(6, 3);
            check_lu_fact_transposed(3, 6, 4, At, 0);
        }

    } // namespace fatrop::test
} // namespace fatrop
