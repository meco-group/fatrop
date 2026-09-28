// Copyright (c) 2026 Lander Vanroye, KU Leuven
// Matrix-RHS overloads reuse the native augmented-system factorization.
#include "fatrop/dense/aug_system_solver.hpp"
#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/linear_algebra/batched_rhs.hpp"

using namespace fatrop;

void AugSystemSolver<DenseType>::allocate_rhs_workspace()
{
    rhs_many_work_.reserve(7);
    rhs_many_work_.emplace_back(v_RSQrqt_tilde_.m(), rhs_capacity_);
    rhs_many_work_.emplace_back(v_Ggt_stripe_.m(), rhs_capacity_);
    rhs_many_work_.emplace_back(v_Ggt_tilde_.m(), rhs_capacity_);
    rhs_many_work_.emplace_back(v_GgLt_.m(), rhs_capacity_);
    rhs_many_work_.emplace_back(v_RSQrqt_hat_.m(), rhs_capacity_);
    rhs_many_work_.emplace_back(rhs_capacity_, v_Llt_.m());
    rhs_many_work_.emplace_back(v_Ggt_ineq_.m(), rhs_capacity_);
}

namespace
{
// BLASFEO optimizes right TRSM for skinny batches. Transpose once around both
// triangular solves, retaining a column-per-RHS API and a cached workspace.
void solve_cholesky_negative(Index n, const MatRealView &L, const MatRealView &R, Index ri,
                             MatRealView &X, Index xi, MatRealView &scratch)
{
    if (n == 0)
        return;
    getr(n, R.n(), R, ri, 0, scratch, 0, 0);
    trsm_rltn(R.n(), n, 1.0, L, 0, 0, scratch, 0, 0, scratch, 0, 0);
    trsm_rlnn(R.n(), n, -1.0, L, 0, 0, scratch, 0, 0, scratch, 0, 0);
    getr(R.n(), n, scratch, 0, 0, X, xi, 0);
}
} // namespace

LinsolReturnFlag AugSystemSolver<DenseType>::solve_rhs(const ProblemInfo<DenseType> &info,
                                                       const Jacobian<DenseType> &jacobian,
                                                       const Hessian<DenseType> &hessian,
                                                       const VecRealView &D_s,
                                                       const MatRealView &f, const MatRealView &g,
                                                       MatRealView &x, MatRealView &eq_mult)
{
    batched::check_dimensions(info.number_of_tangent_variables, info.number_of_eq_constraints,
                              f, g, x, eq_mult);
    if (f.n() == 0) return LinsolReturnFlag::SUCCESS;
    if (f.n() > rhs_capacity_)
        return batched::solve_panels(rhs_capacity_, f, g, x, eq_mult,
            [&](const MatRealView &fp, const MatRealView &gp, MatRealView &xp, MatRealView &lp)
            { return solve_rhs(info, jacobian, hessian, D_s, fp, gp, xp, lp); });
    auto m_RSQrqt_tilde_ = rhs_many_work_[0].block(rhs_many_work_[0].m(), f.n(), 0, 0);
    auto m_Ggt_stripe_ = rhs_many_work_[1].block(rhs_many_work_[1].m(), f.n(), 0, 0);
    auto m_Ggt_tilde_ = rhs_many_work_[2].block(rhs_many_work_[2].m(), f.n(), 0, 0);
    auto m_GgLt_ = rhs_many_work_[3].block(rhs_many_work_[3].m(), f.n(), 0, 0);
    auto m_RSQrqt_hat_ = rhs_many_work_[4].block(rhs_many_work_[4].m(), f.n(), 0, 0);
    auto m_Llt_transposed = rhs_many_work_[5].block(f.n(), rhs_many_work_[5].n(), 0, 0);
    auto m_Ggt_ineq_ = rhs_many_work_[6].block(rhs_many_work_[6].m(), f.n(), 0, 0);

    const Index nx = info.dims.nx_tangent;
    const Index ng = info.dims.ng;
    const Index ng_ineq = info.dims.ng_ineq;

    // Eliminate inequality multipliers from every RHS.
    batched::gecp(nx, f, 0, m_RSQrqt_tilde_, 0);
    if (ng_ineq > 0)
    {
        batched::divide_rows(ng_ineq, g, info.offset_g_eq_slack, D_s, 0, m_Ggt_ineq_, 0);
        batched::gemm_nn(nx, ng_ineq, 1.0, jacobian.Gg_ineqt, 0, 0, m_Ggt_ineq_, 0, 1.0, m_RSQrqt_tilde_, 0,
               m_RSQrqt_tilde_, 0);
    }

    // Reuse the equality-Jacobian LU and reduced-Hessian Cholesky factors.
    if (ng > 0)
    {
        batched::gecp(ng, g, 0, m_Ggt_stripe_, 0);
        batched::rowpe(Pl_, rank_eq_, m_Ggt_stripe_, 0);
        batched::trsm_lutu(rank_eq_, Ggt_stripe_, 0, 0, m_Ggt_stripe_, 0, m_Ggt_stripe_, 0);

        batched::gecpsc(rank_eq_, -1.0, m_Ggt_stripe_, 0, m_Ggt_tilde_, 0);
        batched::trsm_lltn(rank_eq_, Ggt_stripe_, 0, 0, m_Ggt_tilde_, 0, m_Ggt_tilde_, 0);

        batched::rowpe(Pr_, rank_eq_, m_RSQrqt_tilde_, 0);

        batched::gecp(nx, m_RSQrqt_tilde_, 0, m_GgLt_, 0);
        batched::gemm_nn(nx, rank_eq_, 1.0, Ppt_, 0, 0, m_Ggt_tilde_, 0, 1.0, m_GgLt_, 0, m_GgLt_, 0);

        batched::gemm_nn(nx - rank_eq_, rank_eq_, 1.0, Ggt_tilde_, 0, 0, m_GgLt_, 0, 1.0, m_GgLt_, rank_eq_,
               m_RSQrqt_hat_, 0);
        solve_cholesky_negative(nx - rank_eq_, Llt_, m_RSQrqt_hat_, 0, x, rank_eq_, m_Llt_transposed);
        batched::gecp(rank_eq_, m_Ggt_tilde_, 0, x, 0);
        batched::gemm_tn(nx - rank_eq_, rank_eq_, 1.0, Ggt_tilde_, 0, 0, x,
               rank_eq_, 1.0, x, 0, x,
               0);
        batched::gecpsc(rank_eq_, -1.0, m_RSQrqt_tilde_, 0, eq_mult, 0);
        batched::gemm_tn(nx, rank_eq_, -1.0, Ppt_, 0, 0, x, 0, 1.0, eq_mult,
               0, eq_mult, 0);
        batched::trsm_llnn(rank_eq_, Ggt_stripe_, 0, 0, eq_mult, 0, eq_mult,
                 0);
        batched::trsm_lunu(rank_eq_, rank_eq_, Ggt_stripe_, 0, 0, eq_mult, 0, eq_mult,
                 0);
        batched::rowpei(Pl_, rank_eq_, eq_mult, 0);
        batched::rowpei(Pr_, rank_eq_, x, 0);
    }
    else
    {
        solve_cholesky_negative(nx, Llt_, m_RSQrqt_tilde_, 0, x, 0, m_Llt_transposed);
    }

    if (ng_ineq > 0)
    {
        batched::gemm_tn(nx, ng_ineq, 1.0, jacobian.Gg_ineqt, 0, 0, x, 0, 1.0, g,
               info.offset_g_eq_slack, eq_mult, info.offset_g_eq_slack);
        batched::divide_rows(ng_ineq, eq_mult, info.offset_g_eq_slack, D_s, 0, eq_mult, info.offset_g_eq_slack);
    }
    return LinsolReturnFlag::SUCCESS;
}

LinsolReturnFlag AugSystemSolver<DenseType>::solve_rhs(const ProblemInfo<DenseType> &info,
                                                       const Jacobian<DenseType> &jacobian,
                                                       const Hessian<DenseType> &hessian,
                                                       const VecRealView &D_eq,
                                                       const VecRealView &D_s,
                                                       const MatRealView &f, const MatRealView &g,
                                                       MatRealView &x, MatRealView &eq_mult)
{
    batched::check_dimensions(info.number_of_tangent_variables, info.number_of_eq_constraints,
                              f, g, x, eq_mult);
    if (f.n() == 0) return LinsolReturnFlag::SUCCESS;
    if (f.n() > rhs_capacity_)
        return batched::solve_panels(rhs_capacity_, f, g, x, eq_mult,
            [&](const MatRealView &fp, const MatRealView &gp, MatRealView &xp, MatRealView &lp)
            { return solve_rhs(info, jacobian, hessian, D_eq, D_s, fp, gp, xp, lp); });
    auto m_RSQrqt_tilde_ = rhs_many_work_[0].block(rhs_many_work_[0].m(), f.n(), 0, 0);
    auto m_Ggt_stripe_ = rhs_many_work_[1].block(rhs_many_work_[1].m(), f.n(), 0, 0);
    auto m_Ggt_tilde_ = rhs_many_work_[2].block(rhs_many_work_[2].m(), f.n(), 0, 0);
    auto m_GgLt_ = rhs_many_work_[3].block(rhs_many_work_[3].m(), f.n(), 0, 0);
    auto m_RSQrqt_hat_ = rhs_many_work_[4].block(rhs_many_work_[4].m(), f.n(), 0, 0);
    auto m_Llt_transposed = rhs_many_work_[5].block(f.n(), rhs_many_work_[5].n(), 0, 0);
    auto m_Ggt_ineq_ = rhs_many_work_[6].block(rhs_many_work_[6].m(), f.n(), 0, 0);

    const Index nx = info.dims.nx_tangent;
    const Index ng = info.dims.ng;
    const Index ng_ineq = info.dims.ng_ineq;

    // Eliminate inequality multipliers from every RHS.
    batched::gecp(nx, f, 0, m_RSQrqt_tilde_, 0);
    if (ng > 0)
    {
        batched::divide_rows(ng, g, 0, D_eq, 0, m_Ggt_stripe_, 0);
        batched::gemm_nn(nx, ng, 1.0, jacobian.Gg_eqt, 0, 0, m_Ggt_stripe_, 0, 1.0, m_RSQrqt_tilde_, 0,
               m_RSQrqt_tilde_, 0);
    }
    if (ng_ineq > 0)
    {
        batched::divide_rows(ng_ineq, g, info.offset_g_eq_slack, D_s, 0, m_Ggt_ineq_, 0);
        batched::gemm_nn(nx, ng_ineq, 1.0, jacobian.Gg_ineqt, 0, 0, m_Ggt_ineq_, 0, 1.0, m_RSQrqt_tilde_, 0,
               m_RSQrqt_tilde_, 0);
    }
    solve_cholesky_negative(nx, Llt_, m_RSQrqt_tilde_, 0, x, 0, m_Llt_transposed);

    if (ng > 0)
    {
        batched::gemm_tn(nx, ng, 1.0, jacobian.Gg_eqt, 0, 0, x, 0, 1.0, g,
               0, eq_mult, 0);
        batched::divide_rows(ng, eq_mult, 0, D_eq, 0, eq_mult, 0);
    }
    if (ng_ineq > 0)
    {
        batched::gemm_tn(nx, ng_ineq, 1.0, jacobian.Gg_ineqt, 0, 0, x, 0, 1.0, g,
               info.offset_g_eq_slack, eq_mult, info.offset_g_eq_slack);
        batched::divide_rows(ng_ineq, eq_mult, info.offset_g_eq_slack, D_s, 0, eq_mult, info.offset_g_eq_slack);
    }
    return LinsolReturnFlag::SUCCESS;
}
