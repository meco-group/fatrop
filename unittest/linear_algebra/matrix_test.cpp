#include "fatrop/context/generic.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include <gtest/gtest.h>

namespace fatrop
{

TEST(MatrixTest, Creation)
{
    EXPECT_NO_THROW({
        MatrixAllocated mat(3, 4);
        EXPECT_EQ(mat.m(), 3);
        EXPECT_EQ(mat.n(), 4);
    });
}

TEST(MatrixTest, ElementAccess)
{
    MatrixAllocated mat(2, 2);
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
    MatrixAllocated mat1(2, 2);
    mat1(0, 0) = 1.0;
    mat1(0, 1) = 2.0;
    mat1(1, 0) = 3.0;
    mat1(1, 1) = 4.0;

    EXPECT_NO_THROW({
        MatrixAllocated mat2(std::move(mat1));
        EXPECT_EQ(mat2.m(), 2);
        EXPECT_EQ(mat2.n(), 2);
        EXPECT_DOUBLE_EQ(mat2(0, 0), 1.0);
        EXPECT_DOUBLE_EQ(mat2(0, 1), 2.0);
        EXPECT_DOUBLE_EQ(mat2(1, 0), 3.0);
        EXPECT_DOUBLE_EQ(mat2(1, 1), 4.0);
    });
}

} // namespace fatrop
