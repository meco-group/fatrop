#include "../allocation_tracker.hpp"
#include "fatrop/dense/aug_system_solver.hpp"
#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/ocp/aug_system_solver.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include "fatrop/graph/aug_system_solver.hpp"
#include "fatrop/graph/hessian.hpp"
#include "fatrop/graph/jacobian.hpp"
#include "fatrop/graph/problem_info.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace fatrop;
namespace alloc = fatrop::test::allocation_tracker;

namespace
{
template <typename Solve> LinsolReturnFlag without_allocations(Solve solve)
{
    alloc::start();
    const auto ret = solve();
    const auto count = alloc::stop();
    EXPECT_EQ(count.calls, 0u);
    EXPECT_EQ(count.bytes, 0u);
    return ret;
}

template <typename P>
void check_native_solves(const ProblemInfo<P> &info, Jacobian<P> &jac, Hessian<P> &hess,
                         Index path_equalities, Index path_offset)
{
    const Index n = info.number_of_tangent_variables, ng = info.number_of_eq_constraints;
    const Index ni = info.number_of_slack_variables, max_width = 25;
    VecRealAllocated dx(n), de(path_equalities), ds(ni), f(n), g(ng), x(n), lambda(ng);
    VecRealAllocated hx(n), atlambda(n), ax(ng);
    MatRealAllocated F(n, max_width), G(ng, max_width), X(n, max_width), L(ng, max_width);
    for (Index j = 0; j < max_width; ++j)
    {
        for (Index i = 0; i < n; ++i) F(i, j) = std::sin(0.1 * i + 0.7 * j);
        for (Index i = 0; i < ng; ++i) G(i, j) = std::cos(0.3 * i - 0.5 * j);
    }
    for (Index i = 0; i < n; ++i) f(i) = F(i, 0);
    for (Index i = 0; i < ng; ++i) g(i) = G(i, 0);
    de = 0.3;
    ds = 0.7;
    for (Index capacity : {3, 12})
    {
        SCOPED_TRACE(capacity);
        AugSystemSolver<P> solver(info, capacity);
        for (bool regularized : {false, true})
        {
            SCOPED_TRACE(regularized);
            for (Index repeat = 0; repeat < 3; ++repeat)
            {
                SCOPED_TRACE(repeat);
                dx = 0.1 * (repeat + 1);
                ASSERT_EQ(without_allocations([&] {
                    return regularized ? solver.solve(info, jac, hess, dx, de, ds, f, g, x, lambda)
                                       : solver.solve(info, jac, hess, dx, ds, f, g, x, lambda);
                }), LinsolReturnFlag::SUCCESS);
                ASSERT_EQ(without_allocations([&] {
                    return regularized ? solver.solve_rhs(info, jac, hess, de, ds, f, g, x, lambda)
                                       : solver.solve_rhs(info, jac, hess, ds, f, g, x, lambda);
                }), LinsolReturnFlag::SUCCESS);
                for (Index width : {0, 1, 3, 7, 12, 25, 1})
                {
                    SCOPED_TRACE(width);
                    auto fp = F.block(n, width, 0, 0), gp = G.block(ng, width, 0, 0);
                    auto xp = X.block(n, width, 0, 0), lp = L.block(ng, width, 0, 0);
                    ASSERT_EQ(without_allocations([&] {
                        return regularized ? solver.solve_rhs(info, jac, hess, de, ds, fp, gp, xp, lp)
                                           : solver.solve_rhs(info, jac, hess, ds, fp, gp, xp, lp);
                    }), LinsolReturnFlag::SUCCESS);
                    ASSERT_EQ(without_allocations([&] {
                        return regularized ? solver.solve(info, jac, hess, dx, de, ds, fp, gp, xp, lp)
                                           : solver.solve(info, jac, hess, dx, ds, fp, gp, xp, lp);
                    }), LinsolReturnFlag::SUCCESS);
                    if (width == 0) continue;
                    // Check the scalar first column and the final batched panel.
                    for (Index col : {0, width - 1})
                    {
                        for (Index i = 0; i < n; ++i) x(i) = X(i, col);
                        for (Index i = 0; i < ng; ++i) lambda(i) = L(i, col);
                        hess.apply_on_right(info, x, 0.0, hx, hx);
                        jac.transpose_apply_on_right(info, lambda, 0.0, atlambda, atlambda);
                        jac.apply_on_right(info, x, 0.0, ax, ax);
                        for (Index i = 0; i < n; ++i)
                            EXPECT_NEAR(hx(i) + dx(i) * x(i) + atlambda(i) + F(i, col), 0., 1e-8);
                        for (Index i = 0; i < ni; ++i)
                            ax(info.offset_g_eq_slack + i) -= ds(i) * lambda(info.offset_g_eq_slack + i);
                        if (regularized)
                            for (Index i = 0; i < path_equalities; ++i)
                                ax(path_offset + i) -= de(i) * lambda(path_offset + i);
                        for (Index i = 0; i < ng; ++i)
                            EXPECT_NEAR(ax(i) + G(i, col), 0., 1e-8);
                    }
                }
            }
        }
    }
}
}

TEST(NativeRhsAllocation, DenseFactorizationAndMatrixRhs)
{
    ProblemDims<DenseType> dims{12, 3, 4};
    ProblemInfo<DenseType> info{dims};
    Jacobian<DenseType> jac{dims};
    Hessian<DenseType> hess{dims};
    hess.Hht = 0.;
    hess.Hht.block(12, 12, 0, 0).diagonal() = 2.;
    jac.Gg_eqt = 0.;
    jac.Gg_eqt.block(3, 3, 0, 0).diagonal() = 1.;
    jac.Gg_ineqt = 0.;
    jac.Gg_ineqt.block(4, 4, 5, 0).diagonal() = 1.;
    check_native_solves(info, jac, hess, 3, 0);
}

TEST(NativeRhsAllocation, OcpFactorizationAndMatrixRhs)
{
    ProblemDims<OcpType> dims{5, {2, 2, 2, 2, 0}, {3, 3, 3, 3, 3},
                              {3, 1, 0, 1, 0}, {1, 2, 0, 1, 2}};
    ProblemInfo<OcpType> info{dims};
    Jacobian<OcpType> jac{dims};
    Hessian<OcpType> hess{dims};
    for (Index k = 0; k < dims.K; ++k)
    {
        const Index nu = dims.number_of_tangent_controls[k], nx = dims.number_of_tangent_states[k];
        const Index nv = nu + nx, ni = dims.number_of_ineq_constraints[k];
        hess.RSQrqt[k] = 0.;
        hess.RSQrqt[k].block(nv, nv, 0, 0).diagonal() = 2.;
        jac.Gg_eqt[k] = 0.;
        jac.Gg_ineqt[k] = 0.;
        for (Index i = 0; i < ni; ++i) jac.Gg_ineqt[k](nu + i, i) = 1.;
        if (k == 0) jac.Gg_eqt[k].block(nx, nx, nu, 0).diagonal() = 1.;
        else if (dims.number_of_eq_constraints[k]) jac.Gg_eqt[k](0, 0) = 1.;
        if (k + 1 < dims.K)
        {
            jac.BAbt[k] = 0.;
            for (Index i = 0; i < nx; ++i) jac.BAbt[k](nu + i, i) = 0.9;
            for (Index i = 0; i < nu; ++i) jac.BAbt[k](i, i) = 0.2;
        }
    }
    check_native_solves(info, jac, hess, info.number_of_g_eq_path, info.offset_g_eq_path);
}

TEST(NativeRhsAllocation, GraphFactorizationAndMatrixRhs)
{
    ProblemDims<GraphProblem> dims{{2, 3, 4, 1}, {{0, 1}, {1, 2}, {2, 3}, {0, 3}}, {1, 2, 0, 1}};
    ProblemInfo<GraphProblem> info{dims};
    Jacobian<GraphProblem> jac{dims};
    Hessian<GraphProblem> hess{dims};
    hess.H.set_zero();
    for (Index k = 0; k < dims.num_blocks; ++k)
    {
        auto h = hess.H.block(k, k);
        h.diagonal() = 3.;
        jac.Gg_ineqt[k] = 0.;
        for (Index i = 0; i < dims.ng_ineq[k]; ++i) jac.Gg_ineqt[k](i, i) = 1.;
    }
    for (const auto &edge : dims.off_diag_edges)
    {
        auto h = hess.H.block(edge.second, edge.first);
        h.diagonal() = -0.1;
    }
    check_native_solves(info, jac, hess, 0, 0);
}
