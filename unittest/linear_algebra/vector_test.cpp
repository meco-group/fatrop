#include "fatrop/linear_algebra/vector.hpp"
#include <cmath>
#include <gtest/gtest.h>

using namespace fatrop;

class VecTest : public ::testing::Test
{
protected:
    VecRealAllocated vec1{5};
    VecRealAllocated vec2{5};

    void SetUp() override
    {
        // Initialize vectors for testing
        for (Index i = 0; i < 5; ++i)
        {
            vec1(i) = i + 1;
            vec2(i) = (i + 1) * 2;
        }
    }
};

TEST_F(VecTest, AccessOperator)
{
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(vec1(i), i + 1);
        EXPECT_EQ(vec2(i), (i + 1) * 2);
    }
}

TEST_F(VecTest, Size)
{
    EXPECT_EQ(vec1.m(), 5);
    EXPECT_EQ(vec2.m(), 5);
}

TEST_F(VecTest, Sum)
{
    EXPECT_EQ(sum(vec1), 15);
    EXPECT_EQ(sum(vec2), 30);
}

TEST_F(VecTest, NormInf)
{
    EXPECT_EQ(norm_inf(vec1), 5);
    EXPECT_EQ(norm_inf(vec2), 10);
}

TEST_F(VecTest, NormL1)
{
    EXPECT_EQ(norm_l1(vec1), 15);
    EXPECT_EQ(norm_l1(vec2), 30);
}

TEST_F(VecTest, NormL2)
{
    EXPECT_NEAR(norm_l2(vec1), std::sqrt(55), 1e-6);
    EXPECT_NEAR(norm_l2(vec2), std::sqrt(220), 1e-6);
}

TEST_F(VecTest, DotProduct) { EXPECT_EQ(dot(vec1, vec2), 110); }

TEST_F(VecTest, VecRealPlusVecReal)
{
    auto result = vec1 + vec2;
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(result(i), (i + 1) * 3);
    }
}

TEST_F(VecTest, VectorSubtraction)
{
    auto result = vec2 - vec1;
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(result(i), i + 1);
    }
}

TEST_F(VecTest, VecRealTimesScalar)
{
    auto result = 2 * vec1;
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(result(i), (i + 1) * 2);
    }
}

TEST_F(VecTest, VecNumericTimesScalarPlusVecNumericTimesScalar)
{
    auto result = (2 * vec1) + (3 * vec2);
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(result(i), 2 * (i + 1) + 3 * ((i + 1) * 2));
    }
}

TEST_F(VecTest, VecNumericPlusVecNumericTimesScalar)
{
    auto result = vec1 + (2 * vec2);
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(result(i), (i + 1) + 2 * ((i + 1) * 2));
    }
}

TEST_F(VecTest, VecRealTimesVecReal)
{
    auto result = vec1 * vec2;
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(result(i), (i + 1) * (i + 1) * 2);
    }
}

TEST_F(VecTest, VectorDivision)
{
    auto result = vec2 / vec1;
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(result(i), 2);
    }
}

TEST_F(VecTest, AbsoluteValue)
{
    VecRealAllocated vec_neg{5};
    for (Index i = 0; i < 5; ++i)
    {
        vec_neg(i) = -1 * (i + 1);
    }
    auto result = abs(vec_neg);
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(result(i), i + 1);
    }
}

TEST_F(VecTest, Logarithm)
{
    auto result = log(vec1);
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_NEAR(result(i), std::log(i + 1), 1e-6);
    }
}

TEST_F(VecTest, Exponential)
{
    auto result = exp(vec1);
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_NEAR(result(i), std::exp(i + 1), 1e-6);
    }
}

TEST_F(VecTest, Sine)
{
    auto result = sin(vec1);
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_NEAR(result(i), std::sin(i + 1), 1e-6);
    }
}

TEST_F(VecTest, Cosine)
{
    auto result = cos(vec1);
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_NEAR(result(i), std::cos(i + 1), 1e-6);
    }
}

TEST_F(VecTest, Block)
{
    auto block = vec1.block(3, 1);
    EXPECT_EQ(block.m(), 3);
    for (Index i = 0; i < 3; ++i)
    {
        EXPECT_EQ(block(i), vec1(i + 1));
    }
}

TEST_F(VecTest, IfElse)
{
    std::vector<bool> condition(5);
    for (Index i = 0; i < 5; ++i)
    {
        condition[i] = (i % 2 == 0) ? 1 : 0;
    }
    auto result = if_else(condition, vec1, vec2);
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(result(i), (i % 2 == 0) ? vec1(i) : vec2(i));
    }
}

TEST_F(VecTest, IfElseLambda)
{
    auto condition = [](Index i) { return i % 2 == 0; };
    auto result = if_else(condition, vec1, vec2);
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(result(i), (i % 2 == 0) ? vec1(i) : vec2(i));
    }
}

TEST_F(VecTest, VecRealScalar)
{
    VecRealScalar scalar_vec(5, 3.0);
    EXPECT_EQ(scalar_vec.m(), 5);
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(scalar_vec(i), 3.0);
    }
}

TEST_F(VecTest, Negation)
{
    auto result = -vec1;
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(result(i), -(i + 1));
    }
}

TEST_F(VecTest, ScalarDivVecReal)
{
    auto result = 10.0 / vec1;
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_NEAR(result(i), 10.0 / (i + 1), 1e-6);
    }
}

TEST_F(VecTest, Assignment)
{
    VecRealAllocated vec3{5};
    vec3 = vec1 + vec2; // Use an expression that returns a VecReal<Derived>
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(vec3(i), vec1(i) + vec2(i));
    }

    vec3 = 3.0;
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(vec3(i), 3.0);
    }

    // Test assignment from different VecReal expressions
    vec3 = 2.0 * vec1;
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(vec3(i), 2.0 * vec1(i));
    }

    vec3 = vec1 * vec2;
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(vec3(i), vec1(i) * vec2(i));
    }
}

TEST_F(VecTest, MoveConstructor)
{
    VecRealAllocated original{5};
    for (Index i = 0; i < 5; ++i)
    {
        original(i) = i + 1;
    }

    VecRealAllocated moved(std::move(original));

    // Check if the moved vector has the correct values
    EXPECT_EQ(moved.m(), 5);
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(moved(i), i + 1);
    }

    // Check if the original vector's internal pointer is nullified
    EXPECT_EQ(original.vec().mem, nullptr);

    // Accessing elements of the original vector should not cause a crash,
    // but the behavior is undefined. We can't test for specific values here.

    // Moving again from the moved vector should work
    VecRealAllocated moved_again(std::move(moved));
    EXPECT_EQ(moved_again.m(), 5);
    for (Index i = 0; i < 5; ++i)
    {
        EXPECT_EQ(moved_again(i), i + 1);
    }
    EXPECT_EQ(moved.vec().mem, nullptr);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
