// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_LINEAR_ALGEBRA_BATCHED_FACTOR_SOLVE_HPP
#define FATROP_LINEAR_ALGEBRA_BATCHED_FACTOR_SOLVE_HPP

#include "fatrop/linear_algebra/batched_rhs.hpp"
#include "fatrop/linear_algebra/vector.hpp"

namespace fatrop::batched
{
// The native factorization paths include a scalar RHS. Use the first matrix
// column for that solve, then apply the resulting factors to the other columns
// together. Storage is owned by the solver and allocated at construction.
struct FactorSolveWorkspace
{
    FactorSolveWorkspace(Index primal, Index dual)
        : f(primal), g(dual), x(primal), lambda(dual) {}

    VecRealAllocated f, g, x, lambda;
};

template <typename FactorFirst, typename SolveRemaining>
LinsolReturnFlag factor_solve(FactorSolveWorkspace &work,
                              const MatRealView &F, const MatRealView &G,
                              MatRealView &X, MatRealView &Lambda,
                              FactorFirst factor_first, SolveRemaining solve_remaining)
{
    check_dimensions(work.f.m(), work.g.m(), F, G, X, Lambda);
    if (F.n() == 0) return LinsolReturnFlag::SUCCESS;

    for (Index i = 0; i < F.m(); ++i) work.f(i) = F(i, 0);
    for (Index i = 0; i < G.m(); ++i) work.g(i) = G(i, 0);
    const auto status = factor_first(work.f, work.g, work.x, work.lambda);
    if (status != LinsolReturnFlag::SUCCESS) return status;
    for (Index i = 0; i < X.m(); ++i) X(i, 0) = work.x(i);
    for (Index i = 0; i < Lambda.m(); ++i) Lambda(i, 0) = work.lambda(i);
    if (F.n() == 1) return LinsolReturnFlag::SUCCESS;

    auto more_f = F.block(F.m(), F.n() - 1, 0, 1);
    auto more_g = G.block(G.m(), G.n() - 1, 0, 1);
    auto more_x = X.block(X.m(), X.n() - 1, 0, 1);
    auto more_lambda = Lambda.block(Lambda.m(), Lambda.n() - 1, 0, 1);
    return solve_remaining(more_f, more_g, more_x, more_lambda);
}
} // namespace fatrop::batched

#endif
