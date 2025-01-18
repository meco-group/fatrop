
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
#include "fatrop/ocp/pd_solver_orig.hpp"
#include "fatrop/ocp/pd_system_orig.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include "fatrop/ocp/type.hpp"
#include <gtest/gtest.h>
#include <vector>
#include "fatrop/common/exception.hpp"

using namespace fatrop;

class PdTest : public ::testing::Test
{
protected:
    // Create OcpDims object
    int K = 10;                                                   // Number of stages
    std::vector<Index> nx = {20, 10, 10, 10, 10, 2, 0, 1, 10, 5}; // State dimensions for each stage
    std::vector<Index> nu = {1, 4, 2, 10, 1, 30, 4, 1, 10, 0};    // Input dimensions for each stage
    // std::vector<Index> ng = std::vector<Index>(K, 0);      // Equality constraints for each stage
    // std::vector<Index> ng_ineq = std::vector<Index>(K, 1); // Inequality constraints for each
    // stage
    std::vector<Index> ng = {9, 3, 4, 3, 4, 0, 1, 0, 1, 5}; // Equality constraints for each

    // std::vector<Index> ng = {9, 3, 4, 3, 4, 0, 1, 0, 1, 5}; // Equality constraints for each
    std::vector<Index> ng_ineq = {0, 5, 10, 0,  0,
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
    VecRealAllocated sl = VecRealAllocated(info.number_of_slack_variables);
    VecRealAllocated su = VecRealAllocated(info.number_of_slack_variables);
    VecRealAllocated zl = VecRealAllocated(info.number_of_slack_variables);
    VecRealAllocated zu = VecRealAllocated(info.number_of_slack_variables);
    VecRealAllocated rhs_cl = VecRealAllocated(info.number_of_slack_variables);
    VecRealAllocated rhs_cu = VecRealAllocated(info.number_of_slack_variables);
    VecRealAllocated mult = VecRealAllocated(info.number_of_eq_constraints);
    VecRealAllocated rhs_x = VecRealAllocated(info.number_of_primal_variables);
    VecRealAllocated rhs_g = VecRealAllocated(info.number_of_eq_constraints);
    VecRealAllocated rhs_s = VecRealAllocated(info.number_of_slack_variables);
    VecRealAllocated D_x = VecRealAllocated(info.number_of_primal_variables);
    VecRealAllocated D_i = VecRealAllocated(info.number_of_slack_variables);
    VecRealAllocated D_eq = VecRealAllocated(info.number_of_g_eq_path);
    MatRealAllocated full_kkt_matrix =
        MatRealAllocated(info.number_of_primal_variables + info.number_of_eq_constraints,
                         info.number_of_primal_variables + info.number_of_eq_constraints);
    // LinearSystem<PdSystemType<OcpType>> lss{info};
    std::shared_ptr<OcpAugSystemSolver> solver = std::make_shared<OcpAugSystemSolver>(info);
    PdSolverOrig<OcpType> pd_solver = PdSolverOrig<OcpType>(info, solver);
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
            D_eq(i) = 1e-8 * (i + 1);
        }
        for (Index i = 0; i < info.number_of_slack_variables; ++i)
        {
            sl(i) = 1. + 0.1 * i;
            su(i) = 1. + 0.2 * i;
            zl(i) = 1. + 0.3 * i;
            zu(i) = 1. + 0.4 * i;
            rhs_cl(i) = 1. + 0.5 * i;
            rhs_cu(i) = 1. + 0.6 * i;
            D_i(i) =  1.0 + 0.7 * i;
        }
    };
};

TEST_F(PdTest, TestSolve)
{
    LinearSystem<PdSystemType<OcpType>> ls(info, jacobian, hessian, D_x, false, D_eq, D_i, sl, su,
                                           zl, zu, rhs_x, rhs_s, rhs_g, rhs_cl, rhs_cu);
    VecRealAllocated x_full(ls.m());
    VecRealAllocated rhs_save(ls.m());
    VecRealAllocated tmp(ls.m());
    ls.get_rhs(rhs_save);
    LinsolReturnFlag ret = pd_solver.solve(ls, x_full);
    EXPECT_EQ(ret, LinsolReturnFlag::SUCCESS);
    ls.set_rhs(rhs_save);
    ls.apply_on_right(x_full, 1.0, rhs_save, tmp);
    EXPECT_NEAR(norm_inf(tmp), 0, 1e-6);
}

TEST_F(PdTest, TestSolveDegen)
{
    LinearSystem<PdSystemType<OcpType>> ls(info, jacobian, hessian, D_x, true, D_eq, D_i, sl, su,
                                           zl, zu, rhs_x, rhs_s, rhs_g, rhs_cl, rhs_cu);
    VecRealAllocated x_full(ls.m());
    VecRealAllocated rhs_save(ls.m());
    VecRealAllocated tmp(ls.m());
    ls.get_rhs(rhs_save);
    LinsolReturnFlag ret = pd_solver.solve(ls, x_full);
    EXPECT_EQ(ret, LinsolReturnFlag::SUCCESS);
    ls.set_rhs(rhs_save);
    ls.apply_on_right(x_full, 1.0, rhs_save, tmp);
    EXPECT_NEAR(norm_inf(tmp), 0, 1e-6);
}
