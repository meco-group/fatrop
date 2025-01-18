//
// Copyright (c) Lander Vanroye, KU Leuven
//
#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/ocp/dims.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include "fatrop/ocp/type.hpp"
#include <gtest/gtest.h>
#include <vector>

using namespace fatrop;

TEST(JacobianTest, ConstructorTest)
{
    // Create OcpDims object
    int K = 5;                                       // Number of stages
    std::vector<Index> nx = {2, 2, 2, 2, 2};      // State dimensions for each stage
    std::vector<Index> nu = {1, 1, 1, 1, 1};         // Input dimensions for each stage
    std::vector<Index> ng = {1, 0, 0, 0, 2};      // Equality constraints for each stage
    std::vector<Index> ng_ineq = {1, 0, 2, 0, 0}; // Inequality constraints for each stage

    ProblemDims<OcpType> dims(K, nx, nu, ng, ng_ineq);

    // Check if Jacobian object can be constructed without throwing an exception
    EXPECT_NO_THROW({ Jacobian<OcpType> jacobian(dims); });
}

TEST(HessianTest, ConstructorTest)
{
    // Create OcpDims object
    int K = 5;                                       // Number of stages
    std::vector<Index> nx = {2, 2, 2, 2, 2};      // State dimensions for each stage
    std::vector<Index> nu = {1, 1, 1, 1, 1};         // Input dimensions for each stage
    std::vector<Index> ng = {1, 0, 0, 0, 2};      // Equality constraints for each stage
    std::vector<Index> ng_ineq = {1, 0, 2, 0, 0}; // Inequality constraints for each stage

    ProblemDims<OcpType> dims(K, nx, nu, ng, ng_ineq);

    // Check if Hessian object can be constructed without throwing an exception
    EXPECT_NO_THROW({ Hessian<OcpType> hessian(dims); });
}

// TEST(JacobianTest, AssertionViolationTest) {
//     // Create OcpDims object with invalid dimensions
//     int K = 3;  // Number of stages
//     std::vector<Index> nx = {2, 2, 2, 2};  // State dimensions for each stage
//     std::vector<Index> nu = {1, 1, 1};     // Input dimensions for each stage
//     std::vector<Index> ng = {4, 3, 3, 2};  // Equality constraints for each stage (violates
//     assertion) std::vector<Index> ng_ineq = {0, 0, 0, 0};  // Inequality constraints for each
//     stage

//     OcpDims dims(K, nx, nu, ng, ng_ineq);

//     // Check if Jacobian constructor throws an exception due to assertion violation
//     EXPECT_ANY_THROW({
//         Jacobian<OcpType> jacobian(dims);
//     });  // Assuming fatrop_assert throws a std::runtime_error
// }

class JacobianTestOperations : public ::testing::Test
{
protected:
    // Create OcpDims object
    int K = 5;                                    // Number of stages
    std::vector<Index> nx = {2, 2, 2, 2, 2};      // State dimensions for each stage
    std::vector<Index> nu = {1, 1, 4, 1, 0};      // Input dimensions for each stage
    std::vector<Index> ng = {1, 0, 5, 0, 2};      // Equality constraints for each stage
    std::vector<Index> ng_ineq = {1, 0, 3, 0, 0}; // Inequality constraints for each stage

    ProblemDims<OcpType> dims{K, nu, nx, ng, ng_ineq};

    ProblemInfo<OcpType> info{dims};
    // Create Jacobian object
    Jacobian<OcpType> jacobian{dims};
    MatRealAllocated full_matrix =
        MatRealAllocated(info.number_of_eq_constraints, info.number_of_primal_variables);
    VecRealAllocated x = VecRealAllocated(info.number_of_primal_variables);
    VecRealAllocated mult = VecRealAllocated(info.number_of_eq_constraints);
    void SetUp()
    {
        // fill the jacobian with random values
        for (Index k = 0; k < info.dims.K - 1; ++k)
        {
            Index nu = info.dims.number_of_controls[k];
            Index nx = info.dims.number_of_states[k];
            Index nx_next = info.dims.number_of_states[k + 1];
            jacobian.BAbt[k] = 1.0 * (k + 1);
            jacobian.BAbt[k].diagonal() = 2.0 * (k + 1);
            jacobian.Gg_eqt[k] = 3.0 * (k + 1);
            jacobian.Gg_eqt[k].diagonal() = 4.0 * (k + 1);
            jacobian.Gg_ineqt[k] = 5.0 * (k + 1);
            jacobian.Gg_ineqt[k].diagonal() = 6.0 * (k + 1);
        }
        // fill the x vector with random values
        for (Index i = 0; i < info.number_of_primal_variables; ++i)
        {
            x(i) = 1.0 * i;
        }
        // fill the mult vector with random values
        for (Index i = 0; i < info.number_of_eq_constraints; ++i)
        {
            mult(i) = 1.0 * i;
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
            full_matrix.block(nx_next, nu + nx, offs_eq_dyn, offs_ux) =
                transpose(jacobian.BAbt[k].block(nu + nx, nx_next, 0, 0));
            full_matrix.block(nx_next, nx_next, offs_eq_dyn, offs_x_next).diagonal() = -1.0;
        }
        // equality path equality constraints
        for (Index k = 0; k < info.dims.K; ++k)
        {
            Index nu = info.dims.number_of_controls[k];
            Index nx = info.dims.number_of_states[k];
            Index ng = info.dims.number_of_eq_constraints[k];
            Index offset_ux = info.offsets_primal_u[k];
            Index offset_g_eq = info.offsets_g_eq_path[k];
            full_matrix.block(ng, nu + nx, offset_g_eq, offset_ux) =
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
            full_matrix.block(ng_ineq, nu + nx, offset_g_ineq, offset_ux) =
                transpose(jacobian.Gg_ineqt[k].block(nu + nx, ng_ineq, 0, 0));
        }
    };
};

TEST_F(JacobianTestOperations, TestInOut)
{
    // put mult as rhs
    jacobian.set_rhs(info, mult);
    VecRealAllocated out = VecRealAllocated(info.number_of_eq_constraints);
    // get mult back from rhs -> out
    jacobian.get_rhs(info, out);
    for (Index i = 0; i < info.number_of_eq_constraints; ++i)
    {
        EXPECT_NEAR(out(i), mult(i), 1e-10);
    }
}
// check Jacobian if full matrix @ x == jacobian.apply_on_right(x)
TEST_F(JacobianTestOperations, ApplyOnRightTestt)
{
    VecRealAllocated out = VecRealAllocated(info.number_of_eq_constraints);
    jacobian.apply_on_right(info, x, 0.0, out, out);
    VecRealAllocated full_matrix_x(info.number_of_eq_constraints);
    gemv_n(info.number_of_eq_constraints, info.number_of_primal_variables, 1.0, full_matrix, 0, 0,
           x, 0, 0.0, full_matrix_x, 0, full_matrix_x, 0);
    for (Index i = 0; i < info.number_of_eq_constraints; ++i)
    {
        EXPECT_NEAR(out(i), full_matrix_x(i), 1e-10);
    }
}
// check Jacobian if full matrix.T @ mult == jacobian.transpose_apply_on_right(mult)
TEST_F(JacobianTestOperations, TransposeApplyOnRightTest)
{
    VecRealAllocated out = VecRealAllocated(info.number_of_primal_variables);
    jacobian.transpose_apply_on_right(info, mult, 0.0, out, out);
    VecRealAllocated full_matrix_t_mult(info.number_of_primal_variables);
    gemv_t(info.number_of_eq_constraints, info.number_of_primal_variables, 1.0, full_matrix, 0, 0,
           mult, 0, 0.0, full_matrix_t_mult, 0, full_matrix_t_mult, 0);
    for (Index i = 0; i < info.number_of_primal_variables; ++i)
    {
        EXPECT_NEAR(out(i), full_matrix_t_mult(i), 1e-10);
    }
}

class HessianTestOperations : public ::testing::Test
{
protected:
    // Create OcpDims object
    int K = 5;                                    // Number of stages
    std::vector<Index> nx = {2, 2, 2, 2, 2};      // State dimensions for each stage
    std::vector<Index> nu = {1, 1, 4, 1, 0};      // Input dimensions for each stage
    std::vector<Index> ng = {1, 0, 5, 0, 2};      // Equality constraints for each stage
    std::vector<Index> ng_ineq = {1, 0, 3, 0, 0}; // Inequality constraints for each stage

    ProblemDims<OcpType> dims{K, nu, nx, ng, ng_ineq};

    ProblemInfo<OcpType> info{dims};
    // Create Hessian object
    Hessian<OcpType> hessian{dims};
    MatRealAllocated full_matrix =
        MatRealAllocated(info.number_of_primal_variables, info.number_of_primal_variables);
    VecRealAllocated x = VecRealAllocated(info.number_of_primal_variables);
    void SetUp()
    {
        full_matrix = 0.;
        // fill the Hessian with random values
        for (Index k = 0; k < dims.K; ++k)
        {
            hessian.RSQrqt[k].diagonal() = 2.0 * (k + 1);
            hessian.RSQrqt[k](0, 1) = 1e-2;
            hessian.RSQrqt[k](1, 0) = 1e-2;
        }
        // fill the x vector with random values
        for (Index i = 0; i < info.number_of_primal_variables; ++i)
        {
            x(i) = 1.0 * i;
        }
        // populate the full matrix
        for (Index k = 0; k < dims.K; k++)
        {
            Index nu = dims.number_of_controls[k];
            Index nx = dims.number_of_states[k];
            Index offs_ux = info.offsets_primal_u[k];
            full_matrix.block(nu + nx, nu + nx, offs_ux, offs_ux) =
                hessian.RSQrqt[k].block(nu + nx, nu + nx, 0, 0);
        }
    };
};
TEST_F(HessianTestOperations, TestInOut)
{
    // put x as rhs
    hessian.set_rhs(info, x);
    VecRealAllocated out = VecRealAllocated(info.number_of_eq_constraints);
    // get x back from rhs -> out
    hessian.get_rhs(info, out);
    for (Index i = 0; i < info.number_of_primal_variables; ++i)
    {
        EXPECT_NEAR(out(i), x(i), 1e-10);
    }
}
// check Hessian if full matrix @ x == hessian.apply_on_right(x)
TEST_F(HessianTestOperations, ApplyOnRightTest)
{
    VecRealAllocated out = VecRealAllocated(info.number_of_primal_variables);
    hessian.apply_on_right(info, x, 0.0, out, out);
    VecRealAllocated full_matrix_x(info.number_of_primal_variables);
    gemv_n(info.number_of_primal_variables, info.number_of_primal_variables, 1.0, full_matrix, 0, 0,
           x, 0, 0.0, full_matrix_x, 0, full_matrix_x, 0);
    for (Index i = 0; i < info.number_of_primal_variables; ++i)
    {
        EXPECT_NEAR(out(i), full_matrix_x(i), 1e-10);
    }
}
