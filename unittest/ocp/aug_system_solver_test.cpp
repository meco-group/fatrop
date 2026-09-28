#include "../batched_aug_solver_check.hpp"
//
// Copyright (c) Lander Vanroye, KU Leuven
//
#include "../random_matrix.hpp"
#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/ocp/aug_system_solver.hpp"
#include "fatrop/ocp/dims.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include "fatrop/ocp/type.hpp"
#include <gtest/gtest.h>
#include <vector>

using namespace fatrop;

class AugSystemSolverTest : public ::testing::Test
{
protected:
    // Create OcpDims object
    int K = 10;                                                   // Number of stages
    std::vector<Index> nx = {20, 10, 10, 10, 10, 2, 0, 1, 10, 5}; // State dimensions for each stage
    std::vector<Index> nu = {1, 4, 2, 10, 1, 30, 4, 1, 10, 0};    // Input dimensions for each stage
    std::vector<Index> ng = {9, 3, 4, 3, 4, 0, 1, 0, 1, 5}; // Equality constraints for each stage
    std::vector<Index> ng_ineq = {0, 5, 10, 0,   0,
                                  0, 0, 0,  10, 0}; // Inequality constraints for each stage

    ProblemDims<OcpType> dims{K, nu, nx, ng, ng_ineq};

    ProblemInfo<OcpType> info{dims};
    // Create Jacobian object
    Jacobian<OcpType> jacobian{dims};
    MatRealAllocated full_matrix_jacobian =
        MatRealAllocated(info.number_of_eq_constraints, info.number_of_primal_variables);
    Hessian<OcpType> hessian{dims};
    MatRealAllocated full_matrix_hessian =
        MatRealAllocated(info.number_of_primal_variables, info.number_of_primal_variables);
    VecRealAllocated x = VecRealAllocated(info.number_of_primal_variables);
    VecRealAllocated mult = VecRealAllocated(info.number_of_eq_constraints);
    VecRealAllocated rhs_x = VecRealAllocated(info.number_of_primal_variables);
    VecRealAllocated rhs_g = VecRealAllocated(info.number_of_eq_constraints);
    VecRealAllocated D_x = VecRealAllocated(info.number_of_primal_variables);
    VecRealAllocated D_s = VecRealAllocated(info.number_of_slack_variables);
    VecRealAllocated D_eq = VecRealAllocated(info.number_of_g_eq_path);
    MatRealAllocated full_kkt_matrix =
        MatRealAllocated(info.number_of_primal_variables + info.number_of_eq_constraints,
                         info.number_of_primal_variables + info.number_of_eq_constraints);
    AugSystemSolver<OcpType> solver = AugSystemSolver<OcpType>(info);
    void SetUp()
    {
        x = 0;
        full_matrix_jacobian = 0.;
        // fill the jacobian with random values
        for (Index k = 0; k < info.dims.K; ++k)
        {
            Index nu = info.dims.number_of_controls[k];
            Index nx = info.dims.number_of_states[k];
            if (k < info.dims.K - 1)
            {
                Index nx_next = info.dims.number_of_states[k + 1];
                jacobian.BAbt[k].block(nu + nx, nx_next, 0, 0) =
                    ::test::random_matrix(nu + nx, nx_next);
            }
            jacobian.Gg_eqt[k].block(nu + nx, info.dims.number_of_eq_constraints[k], 0, 0) =
                ::test::random_matrix(nu + nx, info.dims.number_of_eq_constraints[k]);
            jacobian.Gg_ineqt[k].block(nu + nx, info.dims.number_of_ineq_constraints[k], 0, 0) =
                ::test::random_matrix(nu + nx, info.dims.number_of_ineq_constraints[k]);
        }
        // fill the Hessian with random values
        for (Index k = 0; k < dims.K; ++k)
        {
            Index nu = info.dims.number_of_controls[k];
            Index nx = info.dims.number_of_states[k];
            hessian.RSQrqt[k].block(nu + nx, nu + nx, 0, 0) = ::test::random_spd_matrix(nu + nx);
        }
        // add dynamics constraints
        for (Index k = 0; k < info.dims.K - 1; ++k)
        {
            Index nu = info.dims.number_of_controls[k];
            Index nx = info.dims.number_of_states[k];
            Index offs_ux = info.offsets_primal_u[k];
            Index offs_x_next = info.offsets_primal_x[k + 1];
            Index nx_next = info.dims.number_of_states[k + 1];
            Index offs_eq_dyn = info.offsets_g_eq_dyn[k];
            full_matrix_jacobian.block(nx_next, nu + nx, offs_eq_dyn, offs_ux) =
                transpose(jacobian.BAbt[k].block(nu + nx, nx_next, 0, 0));
            full_matrix_jacobian.block(nx_next, nx_next, offs_eq_dyn, offs_x_next).diagonal() =
                -1.0;
        }
        // equality path equality constraints
        for (Index k = 0; k < info.dims.K; ++k)
        {
            Index nu = info.dims.number_of_controls[k];
            Index nx = info.dims.number_of_states[k];
            Index ng = info.dims.number_of_eq_constraints[k];
            Index offset_ux = info.offsets_primal_u[k];
            Index offset_g_eq = info.offsets_g_eq_path[k];
            full_matrix_jacobian.block(ng, nu + nx, offset_g_eq, offset_ux) =
                transpose(jacobian.Gg_eqt[k].block(nu + nx, ng, 0, 0));
        }
        // inequality path constraints
        for (Index k = 0; k < info.dims.K; ++k)
        {
            Index nu = info.dims.number_of_controls[k];
            Index nx = info.dims.number_of_states[k];
            Index ng_ineq = info.dims.number_of_ineq_constraints[k];
            Index offset_ux = info.offsets_primal_u[k];
            Index offset_g_ineq = info.offsets_g_eq_slack[k];
            full_matrix_jacobian.block(ng_ineq, nu + nx, offset_g_ineq, offset_ux) =
                transpose(jacobian.Gg_ineqt[k].block(nu + nx, ng_ineq, 0, 0));
        }
        full_matrix_hessian = 0.;
        // populate the full matrix
        for (Index k = 0; k < dims.K; k++)
        {
            Index nu = dims.number_of_controls[k];
            Index nx = dims.number_of_states[k];
            Index offs_ux = info.offsets_primal_u[k];
            full_matrix_hessian.block(nu + nx, nu + nx, offs_ux, offs_ux) =
                hessian.RSQrqt[k].block(nu + nx, nu + nx, 0, 0);
        }
        // set up the full KKT matrix
        full_kkt_matrix.block(info.number_of_primal_variables, info.number_of_primal_variables, 0,
                              0) = full_matrix_hessian;
        full_kkt_matrix.block(info.number_of_primal_variables, info.number_of_eq_constraints, 0,
                              info.number_of_primal_variables) = transpose(full_matrix_jacobian);
        full_kkt_matrix.block(info.number_of_eq_constraints, info.number_of_primal_variables,
                              info.number_of_primal_variables, 0) = full_matrix_jacobian;
        // fill the x vector with random values
        for (Index i = 0; i < info.number_of_primal_variables; ++i)
        {
            rhs_x(i) = 1.0 * i;
            D_x = 1.0 * (i + 0.1);
        }
        // fill the mult vector with random values
        for (Index i = 0; i < info.number_of_eq_constraints; ++i)
        {
            rhs_g(i) = 1.0 * i;
        }

        for (Index i = 0; i < info.number_of_g_eq_path; ++i)
        {
            D_eq(i) = 1.0 * (i + 1);
        }
        for (Index i = 0; i < info.number_of_slack_variables; ++i)
        {
            D_s(i) = 1.0 * (i + 0.1);
        }
    };
};

TEST_F(AugSystemSolverTest, TestSolve)
{
    Index ret = solver.solve(info, jacobian, hessian, D_x, D_s, rhs_x, rhs_g, x, mult);
    EXPECT_EQ(ret, LinsolReturnFlag::SUCCESS);
    VecRealAllocated jac_x(info.number_of_eq_constraints);
    jacobian.apply_on_right(info, x, 0.0, jac_x, jac_x);
    VecRealAllocated rhs_gg(info.number_of_eq_constraints);
    rhs_gg = 0.;
    rhs_gg = rhs_gg + rhs_g + jac_x;
    rhs_gg.block(info.number_of_slack_variables, info.offset_g_eq_slack) =
        rhs_gg.block(info.number_of_slack_variables, info.offset_g_eq_slack) -
        D_s * mult.block(info.number_of_slack_variables, info.offset_g_eq_slack);
    VecRealAllocated grad(info.number_of_primal_variables);
    VecRealAllocated tmp(info.number_of_primal_variables);
    grad = 0;
    hessian.apply_on_right(info, x, 0.0, tmp, tmp);
    grad = grad + tmp + D_x * x;
    jacobian.transpose_apply_on_right(info, mult, 0.0, tmp, tmp);
    grad = grad + tmp;
    grad = grad + rhs_x;
    for (Index i = 0; i < info.number_of_eq_constraints; ++i)
    {
        EXPECT_NEAR(rhs_gg(i), 0, 1e-5);
    }
    for (Index i = 0; i < info.number_of_primal_variables; ++i)
    {
        EXPECT_NEAR(grad(i), 0, 1e-5);
    }
}

TEST_F(AugSystemSolverTest, TestSolveRhs)
{
    Index ret = solver.solve(info, jacobian, hessian, D_x, D_s, rhs_x, rhs_g, x, mult);
    EXPECT_TRUE(ret == LinsolReturnFlag::SUCCESS);
    solver.solve_rhs(info, jacobian, hessian, D_s, rhs_x, rhs_g, x, mult);
    VecRealAllocated jac_x(info.number_of_eq_constraints);
    jacobian.apply_on_right(info, x, 0.0, jac_x, jac_x);
    VecRealAllocated rhs_gg(info.number_of_eq_constraints);
    rhs_gg = 0.;
    rhs_gg = rhs_gg + rhs_g + jac_x;
    rhs_gg.block(info.number_of_slack_variables, info.offset_g_eq_slack) =
        rhs_gg.block(info.number_of_slack_variables, info.offset_g_eq_slack) -
        D_s * mult.block(info.number_of_slack_variables, info.offset_g_eq_slack);
    VecRealAllocated grad(info.number_of_primal_variables);
    VecRealAllocated tmp(info.number_of_primal_variables);
    grad = 0;
    hessian.apply_on_right(info, x, 0.0, tmp, tmp);
    grad = grad + tmp + D_x * x;
    jacobian.transpose_apply_on_right(info, mult, 0.0, tmp, tmp);
    grad = grad + tmp;
    grad = grad + rhs_x;
    for (Index i = 0; i < info.number_of_eq_constraints; ++i)
    {
        EXPECT_NEAR(rhs_gg(i), 0, 1e-5);
    }
    for (Index i = 0; i < info.number_of_primal_variables; ++i)
    {
        EXPECT_NEAR(grad(i), 0, 1e-5);
    }
}

TEST_F(AugSystemSolverTest, TestSolveDegen)
{
    Index ret = solver.solve(info, jacobian, hessian, D_x, D_eq, D_s, rhs_x, rhs_g, x, mult);
    EXPECT_EQ(ret, LinsolReturnFlag::SUCCESS);
    VecRealAllocated jac_x(info.number_of_eq_constraints);
    jacobian.apply_on_right(info, x, 0.0, jac_x, jac_x);
    VecRealAllocated rhs_gg(info.number_of_eq_constraints);
    rhs_gg = 0.;
    rhs_gg = rhs_gg + rhs_g + jac_x;
    rhs_gg.block(info.number_of_slack_variables, info.offset_g_eq_slack) =
        rhs_gg.block(info.number_of_slack_variables, info.offset_g_eq_slack) -
        D_s * mult.block(info.number_of_slack_variables, info.offset_g_eq_slack);
    rhs_gg.block(info.number_of_g_eq_path, info.offset_g_eq_path) =
        rhs_gg.block(info.number_of_g_eq_path, info.offset_g_eq_path) -
        D_eq * mult.block(info.number_of_g_eq_path, info.offset_g_eq_path);
    VecRealAllocated grad(info.number_of_primal_variables);
    VecRealAllocated tmp(info.number_of_primal_variables);
    grad = 0;
    hessian.apply_on_right(info, x, 0.0, tmp, tmp);
    grad = grad + tmp + D_x * x;
    jacobian.transpose_apply_on_right(info, mult, 0.0, tmp, tmp);
    grad = grad + tmp;
    grad = grad + rhs_x;
    for (Index i = 0; i < info.number_of_eq_constraints; ++i)
    {
        EXPECT_NEAR(rhs_gg(i), 0, 1e-5);
    }
    for (Index i = 0; i < info.number_of_primal_variables; ++i)
    {
        EXPECT_NEAR(grad(i), 0, 1e-5);
    }
}

TEST_F(AugSystemSolverTest, TestSolveDegenRhs)
{
    Index ret = solver.solve(info, jacobian, hessian, D_x, D_eq, D_s, rhs_x, rhs_g, x, mult);
    EXPECT_EQ(ret, LinsolReturnFlag::SUCCESS);
    ret = solver.solve_rhs(info, jacobian, hessian, D_eq, D_s, rhs_x, rhs_g, x, mult);
    EXPECT_EQ(ret, LinsolReturnFlag::SUCCESS);
    VecRealAllocated jac_x(info.number_of_eq_constraints);
    jacobian.apply_on_right(info, x, 0.0, jac_x, jac_x);
    VecRealAllocated rhs_gg(info.number_of_eq_constraints);
    rhs_gg = 0.;
    rhs_gg = rhs_gg + rhs_g + jac_x;
    rhs_gg.block(info.number_of_slack_variables, info.offset_g_eq_slack) =
        rhs_gg.block(info.number_of_slack_variables, info.offset_g_eq_slack) -
        D_s * mult.block(info.number_of_slack_variables, info.offset_g_eq_slack);
    rhs_gg.block(info.number_of_g_eq_path, info.offset_g_eq_path) =
        rhs_gg.block(info.number_of_g_eq_path, info.offset_g_eq_path) -
        D_eq * mult.block(info.number_of_g_eq_path, info.offset_g_eq_path);
    VecRealAllocated grad(info.number_of_primal_variables);
    VecRealAllocated tmp(info.number_of_primal_variables);
    grad = 0;
    hessian.apply_on_right(info, x, 0.0, tmp, tmp);
    grad = grad + tmp + D_x * x;
    jacobian.transpose_apply_on_right(info, mult, 0.0, tmp, tmp);
    grad = grad + tmp;
    grad = grad + rhs_x;
    for (Index i = 0; i < info.number_of_eq_constraints; ++i)
    {
        EXPECT_NEAR(rhs_gg(i), 0, 1e-5);
    }
    for (Index i = 0; i < info.number_of_primal_variables; ++i)
    {
        EXPECT_NEAR(grad(i), 0, 1e-5);
    }
}
TEST_F(AugSystemSolverTest, BatchedRhs)
{
    check_batched_augmented_solve(info, jacobian, hessian, solver, D_x, D_eq, D_s,
                                  false, info.offset_g_eq_path);
}

TEST_F(AugSystemSolverTest, BatchedRegularizedRhs)
{
    check_batched_augmented_solve(info, jacobian, hessian, solver, D_x, D_eq, D_s,
                                  true, info.offset_g_eq_path);
}

TEST_F(AugSystemSolverTest, BatchedRhsWithoutIncreasedAccuracy)
{
    solver.set_increased_accuracy(false);
    check_batched_augmented_solve(info, jacobian, hessian, solver, D_x, D_eq, D_s,
                                  false, info.offset_g_eq_path);
}

TEST(OcpBatchedAugSystemSolverTest, LongHorizon)
{
    const Index K = 80;
    std::vector<Index> nx(K, 3), nu(K, 2), ng(K, 0), ni(K, 1);
    nu.back() = 0;
    ng.front() = 3; // Fix x_0; the rank must propagate to the initial-state solve.
    ProblemDims<OcpType> dims{K, nu, nx, ng, ni};
    ProblemInfo<OcpType> info{dims};
    Jacobian<OcpType> jacobian{dims};
    Hessian<OcpType> hessian{dims};
    for (Index k = 0; k < K; ++k)
    {
        const Index nv = nx[k] + nu[k];
        hessian.RSQrqt[k] = 0.0;
        hessian.RSQrqt[k].block(nv, nv, 0, 0).diagonal() = 2.0;
        jacobian.Gg_ineqt[k] = 0.0;
        jacobian.Gg_ineqt[k](nv - 1, 0) = 1.0;
        if (k + 1 < K)
        {
            jacobian.BAbt[k] = 0.0;
            for (Index i = 0; i < 3; ++i) jacobian.BAbt[k](nu[k] + i, i) = 0.9;
            for (Index i = 0; i < 2; ++i) jacobian.BAbt[k](i, i) = 0.2;
        }
    }
    jacobian.Gg_eqt[0] = 0.0;
    jacobian.Gg_eqt[0].block(3, 3, nu[0], 0).diagonal() = 1.0;
    VecRealAllocated dx(info.number_of_tangent_variables), de(info.number_of_g_eq_path);
    VecRealAllocated ds(info.number_of_slack_variables);
    dx = 0.1;
    de = 0.2;
    ds = 0.7;
    AugSystemSolver<OcpType> solver{info};
    check_batched_augmented_solve(info, jacobian, hessian, solver, dx, de, ds,
                                  false, info.offset_g_eq_path);
    check_batched_augmented_solve(info, jacobian, hessian, solver, dx, de, ds,
                                  true, info.offset_g_eq_path);
}

TEST_F(AugSystemSolverTest, BatchedRegularizedRhsWithoutIncreasedAccuracy)
{
    solver.set_increased_accuracy(false);
    check_batched_augmented_solve(info, jacobian, hessian, solver, D_x, D_eq, D_s,
                                  true, info.offset_g_eq_path);
}

TEST(OcpBatchedAugSystemSolverTest, DualRegularizationWithDependentPathEqualities)
{
    // The two initial path equalities have the same Jacobian. Positive D_eq
    // must make the augmented system solvable while leaving dynamics exact.
    ProblemDims<OcpType> dims{3, {1, 1, 0}, {2, 2, 2}, {2, 0, 0}, {1, 0, 1}};
    ProblemInfo<OcpType> info{dims};
    Jacobian<OcpType> jacobian{dims};
    Hessian<OcpType> hessian{dims};
    for (Index k = 0; k < dims.K; ++k)
    {
        const Index nu = dims.number_of_tangent_controls[k];
        const Index nx = dims.number_of_tangent_states[k];
        const Index nv = nu + nx;
        hessian.RSQrqt[k] = 0.0;
        hessian.RSQrqt[k].block(nv, nv, 0, 0).diagonal() = 2.0;
        jacobian.Gg_eqt[k] = 0.0;
        jacobian.Gg_ineqt[k] = 0.0;
        if (dims.number_of_ineq_constraints[k] > 0)
            jacobian.Gg_ineqt[k](nv - 1, 0) = 1.0;
        if (k + 1 < dims.K)
        {
            jacobian.BAbt[k] = 0.0;
            for (Index i = 0; i < nx; ++i) jacobian.BAbt[k](nu + i, i) = 0.8;
            jacobian.BAbt[k](0, 0) = 0.2;
        }
    }
    jacobian.Gg_eqt[0](1, 0) = 1.0;
    jacobian.Gg_eqt[0](1, 1) = 1.0;
    VecRealAllocated dx(info.number_of_tangent_variables), dual_diagonal(info.number_of_eq_constraints);
    dual_diagonal = 0.0;
    auto de = dual_diagonal.block(info.number_of_g_eq_path, info.offset_g_eq_path);
    VecRealAllocated ds(info.number_of_slack_variables);
    dx = 0.1;
    de(0) = 0.2;
    de(1) = 0.7;
    ds = 0.3;
    for (bool increased_accuracy : {false, true})
    {
        SCOPED_TRACE(increased_accuracy);
        AugSystemSolver<OcpType> solver{info};
        solver.set_increased_accuracy(increased_accuracy);
        VecRealAllocated f(info.number_of_tangent_variables), g(info.number_of_eq_constraints);
        VecRealAllocated x(info.number_of_tangent_variables), lambda(info.number_of_eq_constraints);
        f = 0.0;
        g = 0.0;
        ASSERT_EQ(solver.solve(info, jacobian, hessian, dx, ds, f, g, x, lambda),
                  LinsolReturnFlag::NOFULL_RANK);
        check_batched_augmented_solve(info, jacobian, hessian, solver, dx, de, ds,
                                      true, info.offset_g_eq_path);
    }
}
