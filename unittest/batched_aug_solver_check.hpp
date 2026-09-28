#ifndef FATROP_TEST_BATCHED_AUG_SOLVER_CHECK_HPP
#define FATROP_TEST_BATCHED_AUG_SOLVER_CHECK_HPP

#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/nlp/aug_system_solver.hpp"
#include "fatrop/nlp/fwd.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"
#include <gtest/gtest.h>
#include <cmath>

// Independently check every RHS against the original KKT equations and the scalar solve.
template <typename Problem>
void check_batched_augmented_solve(const fatrop::ProblemInfo<Problem> &info,
                                  fatrop::Jacobian<Problem> &jac,
                                  fatrop::Hessian<Problem> &hess,
                                  fatrop::AugSystemSolver<Problem> &solver,
                                  const fatrop::VecRealView &Dx,
                                  const fatrop::VecRealView &De,
                                  const fatrop::VecRealView &Ds,
                                  bool regularized, fatrop::Index equality_offset)
{
    using namespace fatrop;
    const Index n = info.number_of_tangent_variables;
    const Index ng = info.number_of_eq_constraints;
    VecRealAllocated f(n), g(ng), x(n), lam(ng), hx(n), atlam(n), ax(ng);
    f = 0.0;
    g = 0.0;
    const auto flag = regularized ? solver.solve(info, jac, hess, Dx, De, Ds, f, g, x, lam)
                                  : solver.solve(info, jac, hess, Dx, Ds, f, g, x, lam);
    ASSERT_EQ(flag, LinsolReturnFlag::SUCCESS);
    for (Index width : {0, 1, 3, 7, 12, 25, 3})
    {
        SCOPED_TRACE(width);
        MatRealAllocated fs(n + 4, width + 4), gs(ng + 4, width + 4);
        MatRealAllocated xs(n + 4, width + 4), ls(ng + 4, width + 4);
        xs = 101.0;
        ls = 101.0;
        auto F = fs.block(n, width, 1, 2);
        auto G = gs.block(ng, width, 2, 1);
        auto X = xs.block(n, width, 3, 1);
        auto L = ls.block(ng, width, 1, 3);
        for (Index j = 0; j < width; ++j)
        {
            for (Index i = 0; i < n; ++i) F(i, j) = std::sin(0.37 * i + 0.91 * j);
            for (Index i = 0; i < ng; ++i) G(i, j) = std::cos(0.73 * i - 0.51 * j);
        }
        const auto ret = regularized ? solver.solve_rhs(info, jac, hess, De, Ds, F, G, X, L)
                                     : solver.solve_rhs(info, jac, hess, Ds, F, G, X, L);
        ASSERT_EQ(ret, LinsolReturnFlag::SUCCESS);
        for (Index j = 0; j < width; ++j)
        {
            SCOPED_TRACE(j);
            for (Index i = 0; i < n; ++i) f(i) = F(i, j);
            for (Index i = 0; i < ng; ++i) g(i) = G(i, j);
            if (regularized) solver.solve_rhs(info, jac, hess, De, Ds, f, g, x, lam);
            else solver.solve_rhs(info, jac, hess, Ds, f, g, x, lam);
            for (Index i = 0; i < n; ++i) EXPECT_NEAR(X(i, j), x(i), 1e-7);
            for (Index i = 0; i < ng; ++i) EXPECT_NEAR(L(i, j), lam(i), 1e-7);
            for (Index i = 0; i < n; ++i) x(i) = X(i, j);
            for (Index i = 0; i < ng; ++i) lam(i) = L(i, j);
            hess.apply_on_right(info, x, 0.0, hx, hx);
            jac.transpose_apply_on_right(info, lam, 0.0, atlam, atlam);
            jac.apply_on_right(info, x, 0.0, ax, ax);
            for (Index i = 0; i < n; ++i)
                EXPECT_NEAR(hx(i) + Dx(i) * x(i) + atlam(i) + f(i), 0.0, 1e-7);
            for (Index i = 0; i < Ds.m(); ++i)
                ax(info.offset_g_eq_slack + i) -= Ds(i) * lam(info.offset_g_eq_slack + i);
            if (regularized)
                for (Index i = 0; i < De.m(); ++i)
                    ax(equality_offset + i) -= De(i) * lam(equality_offset + i);
            for (Index i = 0; i < ng; ++i) EXPECT_NEAR(ax(i) + g(i), 0.0, 1e-7);
        }
        for (Index j = 0; j < xs.n(); ++j) EXPECT_EQ(xs(0, j), 101.0);
        for (Index j = 0; j < ls.n(); ++j) EXPECT_EQ(ls(0, j), 101.0);
    }
}
#endif
