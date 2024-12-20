#include "fatrop/linear_algebra/linear_algebra.hpp"
#include <gtest/gtest.h>

namespace fatrop::test
{
    using fatrop::MatrixAllocated;
    using fatrop::VecAllocated;

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

        MatrixAllocated matrix(size, size);
        for (Index i = 0; i < size; ++i)
            for (Index j = 0; j < size; ++j)
                matrix(i, j) = i * size + j + 1;

        perm.apply_on_rows(size, &matrix.mat());

        VecAllocated row(size);
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

        MatrixAllocated matrix(size, size);
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

        VecAllocated vector(size);
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

        VecAllocated vector(size);
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

        VecAllocated vector(size);
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

    TEST(LUFactorizationTest, LU_FACT_transposed)
    {
        const Index m = 5;
        const Index n = 5;
        const Index n_max = 5;
        Index rank;
        MatrixAllocated A(m, n);
        MatrixAllocated At(n, m);
        PermutationMatrix Pl(m);
        PermutationMatrix Pr(n);

        // Generate a random matrix A
        for (Index i = 0; i < m; ++i)
        {
            for (Index j = 0; j < n; ++j)
            {
                A(i, j) = static_cast<Scalar>(rand()) / RAND_MAX;
                At(j, i) = A(i, j); // Transpose
            }
        }

        // Compute LU factorization
        LU_FACT_transposed(m, n, n_max, rank, &At.mat(), Pl, Pr);

        // Extract L and U from At
        MatrixAllocated L(m, m);
        MatrixAllocated U(m, n);
        for (Index i = 0; i < m; ++i)
        {
            for (Index j = 0; j < n; ++j)
            {
                if (i > j)
                    L(i, j) = At(j, i);
                else if (i == j)
                    L(i, j) = 1.0;
                else
                    L(i, j) = 0.0;

                if (i <= j)
                    U(i, j) = At(j, i);
                else
                    U(i, j) = 0.0;
            }
        }

        // Check if A = Pl * L * U * Pr^T
        MatrixAllocated temp1(m, n);
        MatrixAllocated temp2(m, n);
        MatrixAllocated result(m, n);

        // Compute L * U
        for (Index i = 0; i < m; ++i)
        {
            for (Index j = 0; j < n; ++j)
            {
                Scalar sum = 0;
                for (Index k = 0; k < m; ++k)
                {
                    sum += L(i, k) * U(k, j);
                }
                temp1(i, j) = sum;
            }
        }

        // Apply left permutation Pr to L * U
        Pr.apply_on_rows(n, &temp1.mat());

        // Apply right permutation Pl to (Pr * L * U)
        Pl.apply_on_cols(m, &temp1.mat());

        // Compare A with Pl * L * U * Pl^T
        for (Index i = 0; i < m; ++i)
        {
            for (Index j = 0; j < n; ++j)
            {
                EXPECT_NEAR(A(i, j), temp1(i, j), 1e-10);
            }
        }

        // Check the rank
        EXPECT_EQ(rank, std::min(m, n));
    }

    TEST(LinearAlgebraTest, fatrop_gead_transposed)
    {
        const Index m = 3;
        const Index n = 4;
        MatrixAllocated A(n, m); // Note: A is transposed
        MatrixAllocated B(m, n);
        Scalar alpha = 2.0;

        // Initialize A and B
        for (Index i = 0; i < n; ++i)
            for (Index j = 0; j < m; ++j)
                A(i, j) = i * m + j + 1;

        for (Index i = 0; i < m; ++i)
            for (Index j = 0; j < n; ++j)
                B(i, j) = i * n + j + 1;

        // Call gead_transposed
        fatrop_gead_transposed(m, n, alpha, &A.mat(), 0, 0, &B.mat(), 0, 0);

        // Check results
        for (Index i = 0; i < m; ++i)
        {
            for (Index j = 0; j < n; ++j)
            {
                Scalar expected = (i * n + j + 1) + alpha * (j * m + i + 1);
                EXPECT_NEAR(expected, B(i, j), 1e-10);
            }
        }
    }

    TEST(LinearAlgebraTest, fatrop_dtrsv_unu)
    {
        const Index m = 3;
        MatrixAllocated A(m, m);
        VecAllocated b(m);
        VecAllocated x(m);

        // Initialize A as upper triangular
        for (Index i = 0; i < m; ++i)
            for (Index j = i; j < m; ++j)
            {
                A(i, j) = i * m + j + 1;
                if (i == j)
                    A(i, i) = 1.;
            }

        // Initialize b
        for (Index i = 0; i < m; ++i)
            b(i) = i + 1;


        // Call dtrsv_unu
        fatrop_dtrsv_unu(m, m, &A.mat(), 0, 0, &b.vec(), 0, &x.vec(), 0);

        // Check if Ax = b
        VecAllocated Ax(m);
        for (Index i = 0; i < m; ++i)
        {
            for (Index j = i; j < m; ++j)
                Ax(i) += A(i, j) * x(j);
        }

        for (Index i = 0; i < m; ++i)
            EXPECT_NEAR(b(i), Ax(i), 1e-10);
    }

    TEST(LinearAlgebraTest, fatrop_dtrsv_utu)
    {
        const Index m = 3;
        MatrixAllocated A(m, m);
        VecAllocated b(m);
        VecAllocated x(m);

        // Initialize A as upper triangular with unit diagonal
        for (Index i = 0; i < m; ++i)
        {
            for (Index j = 0; j < m; ++j)
            {
                if (i < j)
                    A(i, j) = i * m + j + 1;
                else if (i == j)
                    A(i, j) = 1.0;
                else
                    A(i, j) = 0.0;
            }
        }

        // Initialize b
        for (Index i = 0; i < m; ++i)
            b(i) = i + 1;


        // Call dtrsv_utu
        fatrop_dtrsv_utu(m, &A.mat(), 0, 0, &b.vec(), 0, &x.vec(), 0);

        // Check if A^T x = b
        VecAllocated Ax(m);
        for (Index i = 0; i < m; ++i)
        {
            for (Index j = i; j < m; ++j)
                Ax(j) += A(i, j) * x(i);
        }

        for (Index i = 0; i < m; ++i)
            EXPECT_NEAR(b(i), Ax(i), 1e-10);
    }

} // namespace fatrop::test
