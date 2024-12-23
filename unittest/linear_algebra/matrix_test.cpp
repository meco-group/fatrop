#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include <gtest/gtest.h>
#include <sstream>

namespace fatrop
{

    TEST(MatrixTest, Creation)
    {
        EXPECT_NO_THROW({
            MatRealAllocated mat(3, 4);
            EXPECT_EQ(mat.m(), 3);
            EXPECT_EQ(mat.n(), 4);
        });
    }

    TEST(MatrixTest, ElementAccess)
    {
        MatRealAllocated mat(2, 2);
        EXPECT_NO_THROW({
            mat(0, 0) = 1.0;
            mat(0, 1) = 2.0;
            mat(1, 0) = 3.0;
            mat(1, 1) = 4.0;
        });

        EXPECT_DOUBLE_EQ(mat(0, 0), 1.0);
        EXPECT_DOUBLE_EQ(mat(0, 1), 2.0);
        EXPECT_DOUBLE_EQ(mat(1, 0), 3.0);
        EXPECT_DOUBLE_EQ(mat(1, 1), 4.0);
    }

    TEST(MatrixTest, MoveConstructor)
    {
        MatRealAllocated mat1(2, 2);
        mat1(0, 0) = 1.0;
        mat1(0, 1) = 2.0;
        mat1(1, 0) = 3.0;
        mat1(1, 1) = 4.0;

        EXPECT_NO_THROW({
            MatRealAllocated mat2(std::move(mat1));
            EXPECT_EQ(mat2.m(), 2);
            EXPECT_EQ(mat2.n(), 2);
            EXPECT_DOUBLE_EQ(mat2(0, 0), 1.0);
            EXPECT_DOUBLE_EQ(mat2(0, 1), 2.0);
            EXPECT_DOUBLE_EQ(mat2(1, 0), 3.0);
            EXPECT_DOUBLE_EQ(mat2(1, 1), 4.0);
        });
    }

    TEST(Matrix1DViewTest, RowView)
    {
        MatRealAllocated mat(3, 4);
        for (Index i = 0; i < 3; ++i)
            for (Index j = 0; j < 4; ++j)
                mat(i, j) = i * 4 + j;

        auto row1 = mat.row(1);
        EXPECT_EQ(row1.m(), 4);
        for (Index j = 0; j < 4; ++j)
            EXPECT_DOUBLE_EQ(row1(j), mat(1, j));

        // Test assignment
        row1 = VecRealScalar(4, 5.0);
        for (Index j = 0; j < 4; ++j)
            EXPECT_DOUBLE_EQ(mat(1, j), 5.0);

        // Test block view
        auto block = row1.block(2, 1);
        EXPECT_EQ(block.m(), 2);
        EXPECT_DOUBLE_EQ(block(0), 5.0);
        EXPECT_DOUBLE_EQ(block(1), 5.0);
    }

    TEST(Matrix1DViewTest, ColumnView)
    {
        MatRealAllocated mat(3, 4);
        for (Index i = 0; i < 3; ++i)
            for (Index j = 0; j < 4; ++j)
                mat(i, j) = i * 4 + j;

        auto col2 = mat.col(2);
        EXPECT_EQ(col2.m(), 3);
        for (Index i = 0; i < 3; ++i)
            EXPECT_DOUBLE_EQ(col2(i), mat(i, 2));

        // Test assignment
        col2 = VecRealScalar(3, 7.0);
        for (Index i = 0; i < 3; ++i)
            EXPECT_DOUBLE_EQ(mat(i, 2), 7.0);

        // Test block view
        auto block = col2.block(2, 0);
        EXPECT_EQ(block.m(), 2);
        EXPECT_DOUBLE_EQ(block(0), 7.0);
        EXPECT_DOUBLE_EQ(block(1), 7.0);
    }

    TEST(Matrix1DViewTest, DiagonalView)
    {
        MatRealAllocated mat(3, 3);
        for (Index i = 0; i < 3; ++i)
            for (Index j = 0; j < 3; ++j)
                mat(i, j) = i * 3 + j;

        auto diag = mat.diagonal();
        EXPECT_EQ(diag.m(), 3);
        for (Index i = 0; i < 3; ++i)
            EXPECT_DOUBLE_EQ(diag(i), mat(i, i));

        // Test assignment
        diag = VecRealScalar(3, 9.0);
        for (Index i = 0; i < 3; ++i)
            EXPECT_DOUBLE_EQ(mat(i, i), 9.0);

        // Test block view
        auto block = diag.block(2, 1);
        EXPECT_EQ(block.m(), 2);
        EXPECT_DOUBLE_EQ(block(0), 9.0);
        EXPECT_DOUBLE_EQ(block(1), 9.0);
    }

    TEST(Matrix1DViewTest, ViewAssignment)
    {
        MatRealAllocated mat1(3, 4);
        MatRealAllocated mat2(3, 4);
        for (Index i = 0; i < 3; ++i)
            for (Index j = 0; j < 4; ++j)
            {
                mat1(i, j) = i * 4 + j;
                mat2(i, j) = (i * 4 + j) * 2;
            }

        // Test row assignment
        for (Index j = 0; j < 4; ++j)
            mat1.row(1)(j) = mat2.row(2)(j);
        for (Index j = 0; j < 4; ++j)
            EXPECT_DOUBLE_EQ(mat1(1, j), mat2(2, j));

        // Test column assignment
        for (Index i = 0; i < 3; ++i)
            mat1.col(2)(i) = mat2.col(1)(i);
        for (Index i = 0; i < 3; ++i)
            EXPECT_DOUBLE_EQ(mat1(i, 2), mat2(i, 1));

        // Test diagonal assignment
        MatRealAllocated mat3(3, 3);
        for (Index i = 0; i < 3; ++i)
            mat3.diagonal()(i) = mat2.row(1)(i);
        for (Index i = 0; i < 3; ++i)
            EXPECT_DOUBLE_EQ(mat3(i, i), mat2(1, i));
    }

    // TEST(Matrix1DViewTest, OutOfBoundsAccess)
    // {
    //     MatRealAllocated mat(3, 4);

    //     EXPECT_THROW(mat.row(3), std::exception);
    //     EXPECT_THROW(mat.col(4), std::exception);
    //     EXPECT_THROW(mat.row(0)(4), std::exception);
    //     EXPECT_THROW(mat.col(0)(3), std::exception);
    //     EXPECT_THROW(mat.diagonal()(3), std::exception);
    // }

    TEST(Matrix1DViewTest, ViewToViewAssignment)
    {
        MatRealAllocated mat1(3, 4);
        MatRealAllocated mat2(3, 4);

        for (Index i = 0; i < 3; ++i)
            for (Index j = 0; j < 4; ++j)
            {
                mat1(i, j) = i * 4 + j;
                mat2(i, j) = (i * 4 + j) * 2;
            }

        // Row to row assignment
        for (Index j = 0; j < 4; ++j)
            mat1.row(1)(j) = mat2.row(2)(j);
        for (Index j = 0; j < 4; ++j)
            EXPECT_DOUBLE_EQ(mat1(1, j), mat2(2, j));

        // Column to column assignment
        for (Index i = 0; i < 3; ++i)
            mat1.col(2)(i) = mat2.col(1)(i);
        for (Index i = 0; i < 3; ++i)
            EXPECT_DOUBLE_EQ(mat1(i, 2), mat2(i, 1));

        // Diagonal to diagonal assignment
        MatRealAllocated mat3(3, 3);
        for (Index i = 0; i < 3; ++i)
            mat3.diagonal()(i) = mat2.diagonal()(i);
        for (Index i = 0; i < 3; ++i)
            EXPECT_DOUBLE_EQ(mat3(i, i), mat2(i, i));
    }

    TEST(Matrix1DViewTest, DifferentViewTypeAssignment)
    {
        MatRealAllocated mat(3, 3);
        for (Index i = 0; i < 3; ++i)
            for (Index j = 0; j < 3; ++j)
                mat(i, j) = i * 3 + j;

        // Row to column assignment
        for (Index i = 0; i < 3; ++i)
            mat.col(1)(i) = mat.row(0)(i);
        for (Index i = 0; i < 3; ++i)
            EXPECT_DOUBLE_EQ(mat(i, 1), mat(0, i));

        // Column to diagonal assignment
        for (Index i = 0; i < 3; ++i)
            mat.diagonal()(i) = mat.col(2)(i);
        for (Index i = 0; i < 3; ++i)
            EXPECT_DOUBLE_EQ(mat(i, i), mat(i, 2));

        // Diagonal to row assignment
        for (Index j = 0; j < 3; ++j)
            mat.row(1)(j) = mat.diagonal()(j);
        for (Index j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(mat(1, j), mat(j, j));
    }

    TEST(MatrixNumericTest, Creation)
    {
        EXPECT_NO_THROW({
            MatRealAllocated mat(3, 4);
            EXPECT_EQ(mat.m(), 3);
            EXPECT_EQ(mat.n(), 4);
        });
    }

    TEST(MatrixNumericTest, ElementAccess)
    {
        MatRealAllocated mat(2, 2);
        EXPECT_NO_THROW({
            mat(0, 0) = 1.0;
            mat(0, 1) = 2.0;
            mat(1, 0) = 3.0;
            mat(1, 1) = 4.0;
        });

        EXPECT_DOUBLE_EQ(mat(0, 0), 1.0);
        EXPECT_DOUBLE_EQ(mat(0, 1), 2.0);
        EXPECT_DOUBLE_EQ(mat(1, 0), 3.0);
        EXPECT_DOUBLE_EQ(mat(1, 1), 4.0);
    }

    TEST(MatrixNumericTest, BlockFunctionality)
    {
        MatRealAllocated mat(4, 4);
        for (Index i = 0; i < 4; ++i)
            for (Index j = 0; j < 4; ++j)
                mat(i, j) = i * 4 + j;

        auto block = mat.block(2, 2, 1, 1);
        EXPECT_EQ(block.m(), 2);
        EXPECT_EQ(block.n(), 2);

        for (Index i = 0; i < 2; ++i)
            for (Index j = 0; j < 2; ++j)
                EXPECT_DOUBLE_EQ(block(i, j), mat(i + 1, j + 1));

        // Modify the block and check if the original matrix is updated
        block(0, 0) = 100.0;
        EXPECT_DOUBLE_EQ(mat(1, 1), 100.0);
    }

    TEST(MatrixNumericTest, BlockAssignment)
    {
        MatRealAllocated mat1(4, 4);
        MatRealAllocated mat2(2, 2);

        for (Index i = 0; i < 4; ++i)
            for (Index j = 0; j < 4; ++j)
                mat1(i, j) = i * 4 + j;

        for (Index i = 0; i < 2; ++i)
            for (Index j = 0; j < 2; ++j)
                mat2(i, j) = (i + j) * 10;

        auto block = mat1.block(2, 2, 1, 1);
        block = mat2;

        for (Index i = 0; i < 2; ++i)
            for (Index j = 0; j < 2; ++j)
                EXPECT_DOUBLE_EQ(mat1(i + 1, j + 1), mat2(i, j));
    }

    TEST(MatrixNumericTest, ScalarAssignment)
    {
        MatRealAllocated mat(3, 3);
        mat = 5.0;

        for (Index i = 0; i < 3; ++i)
            for (Index j = 0; j < 3; ++j)
                EXPECT_DOUBLE_EQ(mat(i, j), 5.0);

        auto block = mat.block(2, 2, 0, 0);
        block = 10.0;

        for (Index i = 0; i < 2; ++i)
            for (Index j = 0; j < 2; ++j)
                EXPECT_DOUBLE_EQ(mat(i, j), 10.0);

        for (Index i = 2; i < 3; ++i)
            for (Index j = 0; j < 3; ++j)
                EXPECT_DOUBLE_EQ(mat(i, j), 5.0);

        for (Index j = 2; j < 3; ++j)
            EXPECT_DOUBLE_EQ(mat(0, j), 5.0);
        EXPECT_DOUBLE_EQ(mat(1, 2), 5.0);
    }

} // namespace fatrop
