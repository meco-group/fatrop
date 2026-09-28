// Copyright (c) 2026 Lander Vanroye, KU Leuven
#include "fatrop/graph/aug_system_solver.hpp"
#include "fatrop/graph/hessian.hpp"
#include "fatrop/graph/jacobian.hpp"
#include "fatrop/graph/problem_info.hpp"

using namespace fatrop;

LinsolReturnFlag AugSystemSolver<GraphProblem>::solve(
    const ProblemInfo<GraphProblem> &info, Jacobian<GraphProblem> &jacobian,
    Hessian<GraphProblem> &hessian, const VecRealView &D_x, const VecRealView &D_s,
    const MatRealView &F, const MatRealView &G, MatRealView &X, MatRealView &Lambda)
{
    return batched::factor_solve(factor_many_work_, F, G, X, Lambda,
        [&](const VecRealView &f, const VecRealView &g, VecRealView &x, VecRealView &lambda)
        { return solve(info, jacobian, hessian, D_x, D_s, f, g, x, lambda); },
        [&](const MatRealView &f, const MatRealView &g, MatRealView &x, MatRealView &lambda)
        { return solve_rhs(info, jacobian, hessian, D_s, f, g, x, lambda); });
}

LinsolReturnFlag AugSystemSolver<GraphProblem>::solve(
    const ProblemInfo<GraphProblem> &info, Jacobian<GraphProblem> &jacobian,
    Hessian<GraphProblem> &hessian, const VecRealView &D_x, const VecRealView &D_eq,
    const VecRealView &D_s, const MatRealView &F, const MatRealView &G,
    MatRealView &X, MatRealView &Lambda)
{
    return batched::factor_solve(factor_many_work_, F, G, X, Lambda,
        [&](const VecRealView &f, const VecRealView &g, VecRealView &x, VecRealView &lambda)
        { return solve(info, jacobian, hessian, D_x, D_eq, D_s, f, g, x, lambda); },
        [&](const MatRealView &f, const MatRealView &g, MatRealView &x, MatRealView &lambda)
        { return solve_rhs(info, jacobian, hessian, D_eq, D_s, f, g, x, lambda); });
}
