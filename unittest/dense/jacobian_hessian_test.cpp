//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//
#include "fatrop/dense/dims.hpp"
#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include <gtest/gtest.h>

using namespace fatrop;

TEST(DenseJacobianTest, ConstructorTest)
{
    ProblemDims<DenseType> dims(8, 3, 4);
    EXPECT_NO_THROW({ Jacobian<DenseType> jacobian(dims); });
}

TEST(DenseHessianTest, ConstructorTest)
{
    ProblemDims<DenseType> dims(8, 3, 4);
    EXPECT_NO_THROW({ Hessian<DenseType> hessian(dims); });
}

class DenseJacobianTestOperations : public ::testing::Test
{
protected:
    static constexpr Index nx_ = 8;
    static constexpr Index ng_ = 3;
    static constexpr Index ng_ineq_ = 4;

    ProblemDims<DenseType> dims{nx_, ng_, ng_ineq_};
    ProblemInfo<DenseType> info{dims};
    Jacobian<DenseType> jacobian{dims};
    MatRealAllocated full_matrix =
        MatRealAllocated(info.number_of_eq_constraints, info.number_of_primal_variables);
    VecRealAllocated x = VecRealAllocated(info.number_of_primal_variables);
    VecRealAllocated mult = VecRealAllocated(info.number_of_eq_constraints);

    void SetUp() override
    {
        // Random-ish, deterministic fill.
        for (Index i = 0; i < nx_; ++i)
            for (Index j = 0; j < ng_; ++j)
                jacobian.Gg_eqt(i, j) = 1.0 + 0.1 * (i + 2 * j);
        for (Index i = 0; i < nx_; ++i)
            for (Index j = 0; j < ng_ineq_; ++j)
                jacobian.Gg_ineqt(i, j) = 2.0 + 0.07 * (i + 3 * j);
        for (Index i = 0; i < info.number_of_primal_variables; ++i)
            x(i) = 1.0 * i;
        for (Index i = 0; i < info.number_of_eq_constraints; ++i)
            mult(i) = 1.0 * i;
        // Reference matrix: equality block on top (rows 0..ng), inequality block below.
        full_matrix = 0.;
        full_matrix.block(ng_, nx_, 0, 0) = transpose(jacobian.Gg_eqt.block(nx_, ng_, 0, 0));
        full_matrix.block(ng_ineq_, nx_, info.offset_g_eq_slack, 0) =
            transpose(jacobian.Gg_ineqt.block(nx_, ng_ineq_, 0, 0));
    }
};

TEST_F(DenseJacobianTestOperations, TestInOut)
{
    jacobian.set_rhs(info, mult);
    VecRealAllocated out(info.number_of_eq_constraints);
    jacobian.get_rhs(info, out);
    for (Index i = 0; i < info.number_of_eq_constraints; ++i)
        EXPECT_NEAR(out(i), mult(i), 1e-10);
}

TEST_F(DenseJacobianTestOperations, ApplyOnRight)
{
    VecRealAllocated out(info.number_of_eq_constraints);
    jacobian.apply_on_right(info, x, 0.0, out, out);
    VecRealAllocated ref(info.number_of_eq_constraints);
    gemv_n(info.number_of_eq_constraints, info.number_of_primal_variables, 1.0, full_matrix, 0, 0,
           x, 0, 0.0, ref, 0, ref, 0);
    for (Index i = 0; i < info.number_of_eq_constraints; ++i)
        EXPECT_NEAR(out(i), ref(i), 1e-10);
}

TEST_F(DenseJacobianTestOperations, TransposeApplyOnRight)
{
    VecRealAllocated out(info.number_of_primal_variables);
    jacobian.transpose_apply_on_right(info, mult, 0.0, out, out);
    VecRealAllocated ref(info.number_of_primal_variables);
    gemv_t(info.number_of_eq_constraints, info.number_of_primal_variables, 1.0, full_matrix, 0, 0,
           mult, 0, 0.0, ref, 0, ref, 0);
    for (Index i = 0; i < info.number_of_primal_variables; ++i)
        EXPECT_NEAR(out(i), ref(i), 1e-10);
}

class DenseHessianTestOperations : public ::testing::Test
{
protected:
    static constexpr Index nx_ = 8;
    static constexpr Index ng_ = 3;
    static constexpr Index ng_ineq_ = 4;

    ProblemDims<DenseType> dims{nx_, ng_, ng_ineq_};
    ProblemInfo<DenseType> info{dims};
    Hessian<DenseType> hessian{dims};
    MatRealAllocated full_matrix =
        MatRealAllocated(info.number_of_primal_variables, info.number_of_primal_variables);
    VecRealAllocated x = VecRealAllocated(info.number_of_primal_variables);

    void SetUp() override
    {
        full_matrix = 0.;
        for (Index i = 0; i < nx_; ++i)
        {
            for (Index j = 0; j < nx_; ++j)
                hessian.Hht(i, j) = 0.1 * (i + j);
            hessian.Hht(i, i) = 2.0 + 0.1 * i;
        }
        for (Index i = 0; i < nx_; ++i)
            x(i) = 1.0 * i;
        full_matrix.block(nx_, nx_, 0, 0) = hessian.Hht.block(nx_, nx_, 0, 0);
    }
};

TEST_F(DenseHessianTestOperations, TestInOut)
{
    hessian.set_rhs(info, x);
    VecRealAllocated out(info.number_of_primal_variables);
    hessian.get_rhs(info, out);
    for (Index i = 0; i < info.number_of_primal_variables; ++i)
        EXPECT_NEAR(out(i), x(i), 1e-10);
}

TEST_F(DenseHessianTestOperations, ApplyOnRight)
{
    VecRealAllocated out(info.number_of_primal_variables);
    hessian.apply_on_right(info, x, 0.0, out, out);
    VecRealAllocated ref(info.number_of_primal_variables);
    gemv_n(info.number_of_primal_variables, info.number_of_primal_variables, 1.0, full_matrix, 0,
           0, x, 0, 0.0, ref, 0, ref, 0);
    for (Index i = 0; i < info.number_of_primal_variables; ++i)
        EXPECT_NEAR(out(i), ref(i), 1e-10);
}
