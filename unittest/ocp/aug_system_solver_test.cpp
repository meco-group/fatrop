//
// Copyright (c) Lander Vanroye, KU Leuven
//
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
    int K = 5;                                    // Number of stages
    std::vector<Index> nx = {20, 2, 2, 4, 2};      // State dimensions for each stage
    std::vector<Index> nu = {1, 4, 4, 10, 0};     // Input dimensions for each stage
    std::vector<Index> ng = {21, 0, 0, 0, 0};      // Equality constraints for each stage
    std::vector<Index> ng_ineq = {0, 0, 0, 0, 0}; // Inequality constraints for each stage

    OcpDims dims{K, nu, nx, ng, ng_ineq};

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
    MatRealAllocated full_kkt_matrix =
        MatRealAllocated(info.number_of_primal_variables + info.number_of_eq_constraints,
                         info.number_of_primal_variables + info.number_of_eq_constraints);
    OcpAugSystemSolver solver = OcpAugSystemSolver(info);
    void SetUp()
    {
        x = 0;
        full_matrix_jacobian = 0.;
        // fill the jacobian with random values
        for (Index k = 0; k < info.dims.K; ++k)
        {
            Index nu = info.dims.number_of_controls[k];
            Index nx = info.dims.number_of_states[k];
            Index nx_next = info.dims.number_of_states[k + 1];
            if (k < info.dims.K - 1)
            {
                jacobian.BAbt[k] = 1.0 * (k + 1);
                jacobian.BAbt[k].diagonal() = 2.0 * (k + 1);
            }
            jacobian.Gg_eqt[k] = 0.1 * (k + 1);
            jacobian.Gg_eqt[k].diagonal() = 4.0 * (k + 1);
            jacobian.Gg_ineqt[k] = .02 * (k + 1);
            jacobian.Gg_ineqt[k].diagonal() = 6.0 * (k + 1);
        }
        // fill the Hessian with random values
        for (Index k = 0; k < dims.K; ++k)
        {
            hessian.RSQrqt[k].diagonal() = 2.0 * (k + 1);
            hessian.RSQrqt[k](0, 1) = 1e-2;
            hessian.RSQrqt[k](1, 0) = 1e-2;
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
        }
        // fill the mult vector with random values
        for (Index i = 0; i < info.number_of_eq_constraints; ++i)
        {
            rhs_g(i) = 1.0 * i;
        }
        hessian.set_rhs(info, rhs_x);
        jacobian.set_rhs(info, rhs_g);
    };
};

TEST_F(AugSystemSolverTest, TestSolve)
{
    VecRealAllocated D_s(info.number_of_slack_variables);
    D_s = 10.;
    VecRealAllocated rhs_s(info.number_of_slack_variables);
    rhs_s = -1.;
    Index ret = solver.solve(info, jacobian, hessian, D_s, rhs_s, x, mult);
    // solver.solve_rhs(info, jacobian, hessian, D_s, rhs_x, rhs_g, rhs_s, x, mult);
    std::cout << "ret: " << ret << std::endl;
    VecRealAllocated jac_x(info.number_of_eq_constraints);
    jacobian.apply_on_right(info, x, jac_x);
    // std::cout << "full matrix jacobian: " << full_matrix_jacobian << std::endl;
    // std::cout << "x: " << x << std::endl;
    // std::cout << "mult: " << mult << std::endl;
    // std::cout << "rhs_g: " << rhs_g << std::endl;
    // std::cout << "jac_x: " << jac_x << std::endl;
    std::cout << "rhs_g + jac_x: " << rhs_g + jac_x << std::endl;
    VecRealAllocated grad(info.number_of_primal_variables);
    VecRealAllocated tmp(info.number_of_primal_variables);
    grad = 0;
    hessian.apply_on_right(info, x, tmp);
    grad = grad + tmp;
    jacobian.transpose_apply_on_right(info, mult, tmp);
    grad = grad + tmp;
    std::cout << "grad " << grad + rhs_x << std::endl;
}
