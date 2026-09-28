// Copyright (c) 2026 Lander Vanroye, KU Leuven
// Matrix-RHS overloads reuse the native augmented-system factorization.
#include "fatrop/graph/aug_system_solver.hpp"
#include "fatrop/graph/hessian.hpp"
#include "fatrop/graph/jacobian.hpp"
#include "fatrop/graph/problem_info.hpp"
#include "fatrop/linear_algebra/batched_rhs.hpp"

using namespace fatrop;

LinsolReturnFlag AugSystemSolver<GraphProblem>::solve_rhs(
    const ProblemInfo<GraphProblem> &info, const Jacobian<GraphProblem> &jacobian,
    const Hessian<GraphProblem> &hessian, const VecRealView &D_s, const MatRealView &F,
    const MatRealView &G, MatRealView &X, MatRealView &Lambda)
{
    batched::check_dimensions(info.dims.nx_tan, info.number_of_eq_constraints, F, G, X, Lambda);
    if (F.n() == 0) return LinsolReturnFlag::SUCCESS;
    if (F.n() > rhs_capacity_)
        return batched::solve_panels(rhs_capacity_, F, G, X, Lambda,
            [&](const MatRealView &fp, const MatRealView &gp, MatRealView &xp, MatRealView &lp)
            { return solve_rhs(info, jacobian, hessian, D_s, fp, gp, xp, lp); });
    auto rhs = rhs_many_->block(rhs_many_->m(), F.n(), 0, 0);
    auto scaled = scaled_many_->block(scaled_many_->m(), F.n(), 0, 0);
    batched::gecp(F.m(), F, 0, rhs, 0);
    for (Index k = 0; k < info.dims.num_blocks; ++k)
    {
        const Index nk = info.dims.block_sizes_tan[k];
        const Index ng = info.dims.ng_ineq[k];
        if (ng == 0 || nk == 0) continue;
        const Index off = info.dims.block_offsets_tan[k];
        const Index goff = info.offset_g_eq_slack + info.dims.ng_ineq_offsets[k];
        batched::divide_rows(ng, G, goff, D_s, goff, scaled, 0);
        batched::gemm_nn(nk, ng, 1.0, jacobian.Gg_ineqt[k], 0, 0,
                           scaled, 0, 1.0, rhs, off, rhs, off);
    }
    chol_->solve_rhs(rhs, X);
    for (Index k = 0; k < info.dims.num_blocks; ++k)
    {
        const Index nk = info.dims.block_sizes_tan[k];
        const Index ng = info.dims.ng_ineq[k];
        if (ng == 0) continue;
        const Index off = info.dims.block_offsets_tan[k];
        const Index goff = info.offset_g_eq_slack + info.dims.ng_ineq_offsets[k];
        batched::gemm_tn(nk, ng, 1.0, jacobian.Gg_ineqt[k], 0, 0,
                                     X, off, 1.0, G, goff, Lambda, goff);
        batched::divide_rows(ng, Lambda, goff, D_s, goff, Lambda, goff);
    }
    return LinsolReturnFlag::SUCCESS;
}

LinsolReturnFlag AugSystemSolver<GraphProblem>::solve_rhs(
    const ProblemInfo<GraphProblem> &info, const Jacobian<GraphProblem> &jacobian,
    const Hessian<GraphProblem> &hessian, const VecRealView &, const VecRealView &D_s,
    const MatRealView &F, const MatRealView &G, MatRealView &X, MatRealView &Lambda)
{
    return solve_rhs(info, jacobian, hessian, D_s, F, G, X, Lambda);
}
