//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/aug_system_solver.hpp"

#include "fatrop/graph/hessian.hpp"
#include "fatrop/graph/jacobian.hpp"
#include "fatrop/graph/problem_info.hpp"
#include "fatrop/linear_algebra/blasfeo_operations.hpp"

using namespace fatrop;

AugSystemSolver<GraphProblem>::AugSystemSolver(const ProblemInfo<GraphProblem> &info)
    : info_(info), M_(*info.sparsity), rhs_(info.dims.nx_tan),
      chol_(std::make_unique<BlockCholeskySolver>(*info.sparsity))
{
}

void AugSystemSolver<GraphProblem>::register_options(OptionRegistry & /*registry*/)
{
    // No tunables exposed at the moment; pivot tolerance is set on the
    // BlockCholeskySolver directly.
}

void AugSystemSolver<GraphProblem>::set_pivot_tol(const Scalar &value)
{
    chol_->set_pivot_tol(value);
}

namespace
{
// Assemble M = H + sum_k A_k^T D_s_k^{-1} A_k + diag(D_x) by adding
// contributions onto a freshly-copied lower-triangle of H.
void assemble_M(const ProblemInfo<GraphProblem> &info, Hessian<GraphProblem> &hessian,
                Jacobian<GraphProblem> &jacobian, const VecRealView &D_x, const VecRealView &D_s,
                BlockPdMatrix &M)
{
    const auto &sp = M.sparsity();
    // 1. Copy the lower-triangle blocks of H into M.
    for (Index j = 0; j < sp.num_blocks(); ++j)
    {
        for (Index i : sp.column_pattern(j))
        {
            MatRealView dst = M.block(i, j);
            MatRealView src = hessian.H.block(i, j);
            gecp(dst.m(), dst.n(), src, 0, 0, dst, 0, 0);
        }
    }

    // 2. Add block-local penalty A_k^T D_s_k^{-1} A_k to the (k,k) diagonal
    //    block. The transposed Jacobian is stored as (n_k + 1) x ng_ineq[k];
    //    only the top n_k rows carry the Jacobian columns.
    for (Index k = 0; k < info.dims.num_blocks; ++k)
    {
        const Index nk = info.dims.block_sizes_tan[k];
        const Index ngk = info.dims.ng_ineq[k];
        if (ngk == 0 || nk == 0)
            continue;
        const Index s_off = info.offset_g_eq_slack + info.dims.ng_ineq_offsets[k];
        // Use the Jacobian's own buffer as a scratch — we scale a copy in
        // place; refactorisations rebuild from scratch each solve so it is
        // safe to mutate.
        MatRealAllocated &Ak = jacobian.Gg_ineqt[k];
        // Build temporary scaled copy.
        MatRealAllocated scaled(nk, ngk);
        gese(nk, ngk, 0.0, scaled, 0, 0);
        // scaled = Ak^T (top nk rows of the transposed-Jacobian matrix).
        gecp(nk, ngk, Ak, 0, 0, scaled, 0, 0);
        for (Index i = 0; i < ngk; ++i)
        {
            const Scalar di = D_s(s_off + i);
            colsc(nk, 1.0 / di, scaled, 0, i);
        }
        // Diagonal block update: M_kk += scaled * Ak^T (lower triangle).
        MatRealView Mkk = M.block(k, k);
        syrk_ln_mn(nk, nk, ngk, 1.0, scaled, 0, 0, Ak, 0, 0, 1.0, Mkk, 0, 0, Mkk, 0, 0);
    }

    // 3. Add diag(D_x) to each diagonal block (tangent-space damping).
    for (Index k = 0; k < info.dims.num_blocks; ++k)
    {
        const Index nk = info.dims.block_sizes_tan[k];
        if (nk == 0)
            continue;
        const Index off = info.dims.block_offsets_tan[k];
        MatRealView Mkk = M.block(k, k);
        diaad(nk, 1.0, D_x, off, Mkk, 0, 0);
    }
}

// rhs_x = f + sum_k A_k^T D_s_k^{-1} g_k.
void assemble_rhs(const ProblemInfo<GraphProblem> &info, const Jacobian<GraphProblem> &jacobian,
                  const VecRealView &D_s, const VecRealView &f, const VecRealView &g,
                  VecRealView &rhs)
{
    veccp(info.dims.nx_tan, f, 0, rhs, 0);
    for (Index k = 0; k < info.dims.num_blocks; ++k)
    {
        const Index nk = info.dims.block_sizes_tan[k];
        const Index ngk = info.dims.ng_ineq[k];
        if (ngk == 0 || nk == 0)
            continue;
        const Index x_off = info.dims.block_offsets_tan[k];
        const Index s_off = info.offset_g_eq_slack + info.dims.ng_ineq_offsets[k];
        // Build the per-block scaled RHS vector g_k / D_s_k, then accumulate
        // A_k^T (g_k / D_s_k) into rhs[block k].
        VecRealAllocated scaled(ngk);
        for (Index i = 0; i < ngk; ++i)
            scaled(i) = g(s_off + i) / D_s(s_off + i);
        gemv_n(nk, ngk, 1.0, jacobian.Gg_ineqt[k], 0, 0, scaled, 0, 1.0, rhs, x_off, rhs, x_off);
    }
}
} // namespace

LinsolReturnFlag AugSystemSolver<GraphProblem>::solve(const ProblemInfo<GraphProblem> &info,
                                                       Jacobian<GraphProblem> &jacobian,
                                                       Hessian<GraphProblem> &hessian,
                                                       const VecRealView &D_x,
                                                       const VecRealView &D_s,
                                                       const VecRealView &f, const VecRealView &g,
                                                       VecRealView &x, VecRealView &eq_mult)
{
    M_.set_zero();
    assemble_M(info, hessian, jacobian, D_x, D_s, M_);
    assemble_rhs(info, jacobian, D_s, f, g, rhs_);

    // The block-Cholesky solver expects to solve M x = -rhs. Wrap rhs_ in a
    // linear system view and call solve_in_place: it will overwrite the rhs_
    // vector with the (positive) solution x.
    LinearSystem<GraphType> ls(M_, rhs_);
    // We need solve_once_impl, not solve_in_place (no iterative refinement —
    // the assembled SPD system here is exact; iterative refinement is the
    // PdSolverOrig's job at the next level up).
    LinsolReturnFlag ret = chol_->solve_once_impl(ls, rhs_);
    if (ret != LinsolReturnFlag::SUCCESS)
        return ret;
    // After solve_once_impl, `rhs_` is left untouched (the result is stashed
    // into rhs_view inside the LinearSystem). Copy out to x.
    veccp(info.dims.nx_tan, rhs_, 0, x, 0);

    // Recover inequality multipliers: lambda_i = D_s^{-1} (A_i x + g_i).
    for (Index k = 0; k < info.dims.num_blocks; ++k)
    {
        const Index nk = info.dims.block_sizes_tan[k];
        const Index ngk = info.dims.ng_ineq[k];
        if (ngk == 0)
            continue;
        const Index x_off = info.dims.block_offsets_tan[k];
        const Index s_off = info.offset_g_eq_slack + info.dims.ng_ineq_offsets[k];
        gemv_t(nk, ngk, 1.0, jacobian.Gg_ineqt[k], 0, 0, x, x_off, 1.0, g, s_off, eq_mult, s_off);
        eq_mult.block(ngk, s_off) = eq_mult.block(ngk, s_off) / D_s.block(ngk, s_off);
    }
    return LinsolReturnFlag::SUCCESS;
}

LinsolReturnFlag AugSystemSolver<GraphProblem>::solve(const ProblemInfo<GraphProblem> &info,
                                                       Jacobian<GraphProblem> &jacobian,
                                                       Hessian<GraphProblem> &hessian,
                                                       const VecRealView &D_x,
                                                       const VecRealView & /*D_eq*/,
                                                       const VecRealView &D_s,
                                                       const VecRealView &f, const VecRealView &g,
                                                       VecRealView &x, VecRealView &eq_mult)
{
    // Graph problems have no equality constraints, so D_eq plays no role.
    return solve(info, jacobian, hessian, D_x, D_s, f, g, x, eq_mult);
}

LinsolReturnFlag AugSystemSolver<GraphProblem>::solve_rhs(const ProblemInfo<GraphProblem> &info,
                                                           const Jacobian<GraphProblem> &jacobian,
                                                           const Hessian<GraphProblem> & /*hessian*/,
                                                           const VecRealView &D_s,
                                                           const VecRealView &f,
                                                           const VecRealView &g, VecRealView &x,
                                                           VecRealView &eq_mult)
{
    // Rebuild only the right-hand side; the Cholesky factor produced by the
    // most recent `solve` call is reused.
    assemble_rhs(info, jacobian, D_s, f, g, rhs_);
    LinearSystem<GraphType> ls(M_, rhs_);
    chol_->solve_rhs_impl(ls, rhs_);
    veccp(info.dims.nx_tan, rhs_, 0, x, 0);
    for (Index k = 0; k < info.dims.num_blocks; ++k)
    {
        const Index nk = info.dims.block_sizes_tan[k];
        const Index ngk = info.dims.ng_ineq[k];
        if (ngk == 0)
            continue;
        const Index x_off = info.dims.block_offsets_tan[k];
        const Index s_off = info.offset_g_eq_slack + info.dims.ng_ineq_offsets[k];
        gemv_t(nk, ngk, 1.0, jacobian.Gg_ineqt[k], 0, 0, x, x_off, 1.0, g, s_off, eq_mult, s_off);
        eq_mult.block(ngk, s_off) = eq_mult.block(ngk, s_off) / D_s.block(ngk, s_off);
    }
    return LinsolReturnFlag::SUCCESS;
}

LinsolReturnFlag AugSystemSolver<GraphProblem>::solve_rhs(const ProblemInfo<GraphProblem> &info,
                                                           const Jacobian<GraphProblem> &jacobian,
                                                           const Hessian<GraphProblem> &hessian,
                                                           const VecRealView & /*D_eq*/,
                                                           const VecRealView &D_s,
                                                           const VecRealView &f,
                                                           const VecRealView &g, VecRealView &x,
                                                           VecRealView &eq_mult)
{
    return solve_rhs(info, jacobian, hessian, D_s, f, g, x, eq_mult);
}
