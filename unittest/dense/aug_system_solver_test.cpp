//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//
#include "../random_matrix.hpp"
#include "fatrop/dense/aug_system_solver.hpp"
#include "fatrop/dense/dims.hpp"
#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include <gtest/gtest.h>

using namespace fatrop;

class DenseAugSystemSolverTest : public ::testing::Test
{
protected:
    static constexpr Index nx_ = 12;
    static constexpr Index ng_ = 4;
    static constexpr Index ng_ineq_ = 5;

    ProblemDims<DenseType> dims{nx_, ng_, ng_ineq_};
    ProblemInfo<DenseType> info{dims};
    Jacobian<DenseType> jacobian{dims};
    Hessian<DenseType> hessian{dims};
    MatRealAllocated full_matrix_jacobian{info.number_of_eq_constraints,
                                          info.number_of_primal_variables};
    MatRealAllocated full_matrix_hessian{info.number_of_primal_variables,
                                         info.number_of_primal_variables};
    VecRealAllocated x{info.number_of_primal_variables};
    VecRealAllocated mult{info.number_of_eq_constraints};
    VecRealAllocated rhs_x{info.number_of_primal_variables};
    VecRealAllocated rhs_g{info.number_of_eq_constraints};
    VecRealAllocated D_x{info.number_of_primal_variables};
    VecRealAllocated D_s{info.number_of_slack_variables};
    VecRealAllocated D_eq{info.number_of_g_eq_slack > 0 ? ng_ : ng_};
    AugSystemSolver<DenseType> solver{info};

    void SetUp() override
    {
        x = 0;
        full_matrix_jacobian = 0.;
        jacobian.Gg_eqt.block(nx_, ng_, 0, 0) = ::test::random_matrix(nx_, ng_);
        jacobian.Gg_ineqt.block(nx_, ng_ineq_, 0, 0) = ::test::random_matrix(nx_, ng_ineq_);
        hessian.Hht.block(nx_, nx_, 0, 0) = ::test::random_spd_matrix(nx_);
        // Assemble the full Jacobian: [A_eq; A_ineq].
        full_matrix_jacobian.block(ng_, nx_, 0, 0) =
            transpose(jacobian.Gg_eqt.block(nx_, ng_, 0, 0));
        full_matrix_jacobian.block(ng_ineq_, nx_, info.offset_g_eq_slack, 0) =
            transpose(jacobian.Gg_ineqt.block(nx_, ng_ineq_, 0, 0));
        full_matrix_hessian = 0.;
        full_matrix_hessian.block(nx_, nx_, 0, 0) = hessian.Hht.block(nx_, nx_, 0, 0);
        for (Index i = 0; i < info.number_of_primal_variables; ++i)
        {
            rhs_x(i) = 1.0 * i;
            D_x(i) = 1.0 * (i + 0.1);
        }
        for (Index i = 0; i < info.number_of_eq_constraints; ++i)
            rhs_g(i) = 1.0 * i;
        for (Index i = 0; i < ng_; ++i)
            D_eq(i) = 1.0 * (i + 1);
        for (Index i = 0; i < info.number_of_slack_variables; ++i)
            D_s(i) = 1.0 * (i + 0.1);
    }

    // Verify the KKT residual: A x + diag(0, D_s) mult + rhs_g_with_penalty = 0
    // and H x + D_x x + A^T mult + rhs_x = 0.
    void check_kkt(bool include_D_eq)
    {
        VecRealAllocated jac_x(info.number_of_eq_constraints);
        jacobian.apply_on_right(info, x, 0.0, jac_x, jac_x);
        VecRealAllocated rhs_gg(info.number_of_eq_constraints);
        rhs_gg = 0.;
        rhs_gg = rhs_gg + rhs_g + jac_x;
        rhs_gg.block(info.number_of_slack_variables, info.offset_g_eq_slack) =
            rhs_gg.block(info.number_of_slack_variables, info.offset_g_eq_slack) -
            D_s * mult.block(info.number_of_slack_variables, info.offset_g_eq_slack);
        if (include_D_eq)
            rhs_gg.block(ng_, 0) = rhs_gg.block(ng_, 0) - D_eq * mult.block(ng_, 0);
        VecRealAllocated grad(info.number_of_primal_variables);
        VecRealAllocated tmp(info.number_of_primal_variables);
        grad = 0.;
        hessian.apply_on_right(info, x, 0.0, tmp, tmp);
        grad = grad + tmp + D_x * x;
        jacobian.transpose_apply_on_right(info, mult, 0.0, tmp, tmp);
        grad = grad + tmp + rhs_x;
        for (Index i = 0; i < info.number_of_eq_constraints; ++i)
            EXPECT_NEAR(rhs_gg(i), 0., 1e-7);
        for (Index i = 0; i < info.number_of_primal_variables; ++i)
            EXPECT_NEAR(grad(i), 0., 1e-7);
    }
};

TEST_F(DenseAugSystemSolverTest, TestSolve)
{
    LinsolReturnFlag ret = solver.solve(info, jacobian, hessian, D_x, D_s, rhs_x, rhs_g, x, mult);
    EXPECT_EQ(ret, LinsolReturnFlag::SUCCESS);
    check_kkt(/*include_D_eq=*/false);
}

TEST_F(DenseAugSystemSolverTest, TestSolveRhs)
{
    LinsolReturnFlag ret = solver.solve(info, jacobian, hessian, D_x, D_s, rhs_x, rhs_g, x, mult);
    ASSERT_EQ(ret, LinsolReturnFlag::SUCCESS);
    ret = solver.solve_rhs(info, jacobian, hessian, D_s, rhs_x, rhs_g, x, mult);
    EXPECT_EQ(ret, LinsolReturnFlag::SUCCESS);
    check_kkt(/*include_D_eq=*/false);
}

TEST_F(DenseAugSystemSolverTest, TestSolveDegen)
{
    LinsolReturnFlag ret =
        solver.solve(info, jacobian, hessian, D_x, D_eq, D_s, rhs_x, rhs_g, x, mult);
    EXPECT_EQ(ret, LinsolReturnFlag::SUCCESS);
    check_kkt(/*include_D_eq=*/true);
}

TEST_F(DenseAugSystemSolverTest, TestSolveDegenRhs)
{
    LinsolReturnFlag ret =
        solver.solve(info, jacobian, hessian, D_x, D_eq, D_s, rhs_x, rhs_g, x, mult);
    ASSERT_EQ(ret, LinsolReturnFlag::SUCCESS);
    ret = solver.solve_rhs(info, jacobian, hessian, D_eq, D_s, rhs_x, rhs_g, x, mult);
    EXPECT_EQ(ret, LinsolReturnFlag::SUCCESS);
    check_kkt(/*include_D_eq=*/true);
}
