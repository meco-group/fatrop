//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//
#include "fatrop/dense/aug_system_solver.hpp"
#include "fatrop/dense/dims.hpp"
#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"

#include <algorithm>

using namespace fatrop;

namespace
{
bool check_reg(const Index m, MAT *sA, const Index ai, const Index aj)
{
    for (Index i = 0; i < m; i++)
    {
        if (blasfeo_matel_wrap(sA, ai + i, aj + i) < 1e-8)
            return false;
    }
    return true;
}
} // namespace

AugSystemSolver<DenseType>::AugSystemSolver(const ProblemInfo<DenseType> &info)
    : Ppt_(info.dims.nx_tangent + 1, info.dims.nx_tangent),
      RSQrqt_tilde_(info.dims.nx_tangent + 1, info.dims.nx_tangent),
      Ggt_stripe_(info.dims.nx_tangent + 1, info.dims.ng),
      Ggt_tilde_(info.dims.nx_tangent + 1, info.dims.nx_tangent),
      GgLt_(info.dims.nx_tangent + 1, info.dims.nx_tangent),
      RSQrqt_hat_(info.dims.nx_tangent + 1, info.dims.nx_tangent),
      Llt_(info.dims.nx_tangent + 1, info.dims.nx_tangent),
      Ggt_ineq_(info.dims.nx_tangent + 1, info.dims.ng_ineq),
      v_RSQrqt_tilde_(info.dims.nx_tangent),
      v_Ggt_stripe_(info.dims.ng),
      v_Ggt_tilde_(info.dims.nx_tangent),
      v_GgLt_(info.dims.nx_tangent),
      v_RSQrqt_hat_(info.dims.nx_tangent),
      v_Llt_(info.dims.nx_tangent),
      v_Ggt_ineq_(info.dims.ng_ineq),
      Pl_(info.dims.ng), Pr_(info.dims.nx_tangent)
{
}

void AugSystemSolver<DenseType>::register_options(OptionRegistry &registry)
{
    registry.register_option("linsol_lu_fact_tol", &AugSystemSolver<DenseType>::set_lu_fact_tol,
                             this);
}

LinsolReturnFlag AugSystemSolver<DenseType>::solve(const ProblemInfo<DenseType> &info,
                                                   Jacobian<DenseType> &jacobian,
                                                   Hessian<DenseType> &hessian,
                                                   const VecRealView &D_x, const VecRealView &D_s,
                                                   const VecRealView &f, const VecRealView &g,
                                                   VecRealView &x, VecRealView &eq_mult)
{
    const Index nx = info.dims.nx_tangent;
    const Index ng = info.dims.ng;
    const Index ng_ineq = info.dims.ng_ineq;

    // ----- 1. Build the augmented Hessian: H + grad row + Σ A_iᵀ Dₛ⁻¹ A_i + diag(D_x).
    rowin(nx, 1.0, f, 0, hessian.Hht, nx, 0);
    gecp(nx + 1, nx, hessian.Hht, 0, 0, RSQrqt_tilde_, 0, 0);

    if (ng_ineq > 0)
    {
        rowin(ng_ineq, 1.0, g, info.offset_g_eq_slack, jacobian.Gg_ineqt, nx, 0);
        gecp(nx + 1, ng_ineq, jacobian.Gg_ineqt, 0, 0, Ggt_ineq_, 0, 0);
        for (Index i = 0; i < ng_ineq; i++)
        {
            const Scalar scaling = 1.0 / D_s(i);
            colsc(nx + 1, scaling, Ggt_ineq_, 0, i);
        }
        syrk_ln_mn(nx + 1, nx, ng_ineq, 1.0, Ggt_ineq_, 0, 0, jacobian.Gg_ineqt, 0, 0, 1.0,
                   RSQrqt_tilde_, 0, 0, RSQrqt_tilde_, 0, 0);
    }
    diaad(nx, 1.0, D_x, 0, RSQrqt_tilde_, 0, 0);
    // Mirror the lower triangle to the upper triangle so that the column
    // permutation below sees a symmetric matrix.
    trtr_l(nx, RSQrqt_tilde_, 0, 0, RSQrqt_tilde_, 0, 0);
    // Copy into Ppt_, which will be progressively reduced.
    gecp(nx + 1, nx, RSQrqt_tilde_, 0, 0, Ppt_, 0, 0);

    // ----- 2. Eliminate equality constraints via a Schur complement.
    if (ng > 0)
    {
        if (ng > nx)
            return LinsolReturnFlag::NOFULL_RANK;

        rowin(ng, 1.0, g, 0, jacobian.Gg_eqt, nx, 0);
        gecp(nx + 1, ng, jacobian.Gg_eqt, 0, 0, Ggt_stripe_, 0, 0);

        // LU factorize the (transposed) equality Jacobian. Pivots are searched in the full
        // nx-column space because there are no controls.
        lu_fact_transposed(ng, nx + 1, nx, rank_eq_, Ggt_stripe_, Pl_, Pr_, lu_fact_tol_);
        if (rank_eq_ < ng)
            return LinsolReturnFlag::NOFULL_RANK;

        // Compute the projection factor Ggt_tilde = -Ggt_stripe[rank:nx+1, :rank] L⁻ᵀ.
        trsm_rlnn(nx - rank_eq_ + 1, rank_eq_, -1.0, Ggt_stripe_, 0, 0, Ggt_stripe_, rank_eq_, 0,
                  Ggt_tilde_, 0, 0);

        // Apply the column permutation to Ppt_ (both rows and columns — symmetric).
        Pr_.apply_on_rows(rank_eq_, &Ppt_.mat());
        Pr_.apply_on_cols(rank_eq_, &Ppt_.mat());

        // GgLt_ = Ppt[rank:nx+1, :nx]  +  Ggt_tilde @ Ppt[:nx, :rank]ᵀ
        gecp(nx - rank_eq_ + 1, nx, Ppt_, rank_eq_, 0, GgLt_, 0, 0);
        gemm_nt(nx - rank_eq_ + 1, nx, rank_eq_, 1.0, Ggt_tilde_, 0, 0, Ppt_, 0, 0, 1.0, GgLt_, 0,
                0, GgLt_, 0, 0);

        // Schur complement: RSQrqt_hat = GgLt[:, :rank] @ Ggt_tildeᵀ + GgLt[:, rank:].
        syrk_ln_mn(nx - rank_eq_ + 1, nx - rank_eq_, rank_eq_, 1.0, GgLt_, 0, 0, Ggt_tilde_, 0, 0,
                   1.0, GgLt_, 0, rank_eq_, RSQrqt_hat_, 0, 0);

        potrf_l_mn(nx - rank_eq_ + 1, nx - rank_eq_, RSQrqt_hat_, 0, 0, Llt_, 0, 0);
        if (!check_reg(nx - rank_eq_, &Llt_.mat(), 0, 0))
            return LinsolReturnFlag::INDEFINITE;
    }
    else
    {
        rank_eq_ = 0;
        potrf_l_mn(nx + 1, nx, Ppt_, 0, 0, Llt_, 0, 0);
        if (!check_reg(nx, &Llt_.mat(), 0, 0))
            return LinsolReturnFlag::INDEFINITE;
    }

    // ----- 3. Forward substitution for the primal direction x.
    if (ng > 0)
    {
        // Free part of x (last nx - rank components).
        rowex(nx - rank_eq_, -1.0, Llt_, nx - rank_eq_, 0, x, rank_eq_);
        trsv_ltn(nx - rank_eq_, Llt_, 0, 0, x, rank_eq_, x,
                 rank_eq_);
        // Constrained part of x (first rank components).
        rowex(rank_eq_, 1.0, Ggt_tilde_, nx - rank_eq_, 0, x, 0);
        gemv_t(nx - rank_eq_, rank_eq_, 1.0, Ggt_tilde_, 0, 0, x,
               rank_eq_, 1.0, x, 0, x,
               0);
        // Equality multipliers.
        rowex(rank_eq_, -1.0, Ppt_, nx, 0, eq_mult, 0);
        gemv_t(nx, rank_eq_, -1.0, Ppt_, 0, 0, x, 0, 1.0, eq_mult,
               0, eq_mult, 0);
        trsv_lnn(rank_eq_, Ggt_stripe_, 0, 0, eq_mult, 0, eq_mult,
                 0);
        trsv_unu(rank_eq_, rank_eq_, Ggt_stripe_, 0, 0, eq_mult, 0, eq_mult,
                 0);
        Pl_.apply_inverse(rank_eq_, &eq_mult.vec(), 0);
        Pr_.apply_inverse(rank_eq_, &x.vec(), 0);
    }
    else
    {
        rowex(nx, -1.0, Llt_, nx, 0, x, 0);
        trsv_ltn(nx, Llt_, 0, 0, x, 0, x, 0);
    }

    // ----- 4. Inequality multipliers: λ_i = D_s⁻¹ (A_i x + g_i).
    if (ng_ineq > 0)
    {
        gemv_t(nx, ng_ineq, 1.0, jacobian.Gg_ineqt, 0, 0, x, 0, 1.0, g,
               info.offset_g_eq_slack, eq_mult, info.offset_g_eq_slack);
        eq_mult.block(ng_ineq, info.offset_g_eq_slack) =
            eq_mult.block(ng_ineq, info.offset_g_eq_slack) / D_s.block(ng_ineq, 0);
    }

    return LinsolReturnFlag::SUCCESS;
}

LinsolReturnFlag AugSystemSolver<DenseType>::solve(const ProblemInfo<DenseType> &info,
                                                   Jacobian<DenseType> &jacobian,
                                                   Hessian<DenseType> &hessian,
                                                   const VecRealView &D_x, const VecRealView &D_eq,
                                                   const VecRealView &D_s, const VecRealView &f,
                                                   const VecRealView &g, VecRealView &x,
                                                   VecRealView &eq_mult)
{
    const Index nx = info.dims.nx_tangent;
    const Index ng = info.dims.ng;
    const Index ng_ineq = info.dims.ng_ineq;

    // ----- 1. Build the regularized augmented Hessian:
    //  H + grad + Σ A_eᵀ D_eq⁻¹ A_e + Σ A_iᵀ D_s⁻¹ A_i + diag(D_x).
    rowin(nx, 1.0, f, 0, hessian.Hht, nx, 0);
    gecp(nx + 1, nx, hessian.Hht, 0, 0, RSQrqt_tilde_, 0, 0);

    if (ng > 0)
    {
        rowin(ng, 1.0, g, 0, jacobian.Gg_eqt, nx, 0);
        gecp(nx + 1, ng, jacobian.Gg_eqt, 0, 0, Ggt_stripe_, 0, 0);
        for (Index i = 0; i < ng; i++)
        {
            const Scalar scaling = 1.0 / D_eq(i);
            colsc(nx + 1, scaling, Ggt_stripe_, 0, i);
        }
        syrk_ln_mn(nx + 1, nx, ng, 1.0, Ggt_stripe_, 0, 0, jacobian.Gg_eqt, 0, 0, 1.0,
                   RSQrqt_tilde_, 0, 0, RSQrqt_tilde_, 0, 0);
    }
    if (ng_ineq > 0)
    {
        rowin(ng_ineq, 1.0, g, info.offset_g_eq_slack, jacobian.Gg_ineqt, nx, 0);
        gecp(nx + 1, ng_ineq, jacobian.Gg_ineqt, 0, 0, Ggt_ineq_, 0, 0);
        for (Index i = 0; i < ng_ineq; i++)
        {
            const Scalar scaling = 1.0 / D_s(i);
            colsc(nx + 1, scaling, Ggt_ineq_, 0, i);
        }
        syrk_ln_mn(nx + 1, nx, ng_ineq, 1.0, Ggt_ineq_, 0, 0, jacobian.Gg_ineqt, 0, 0, 1.0,
                   RSQrqt_tilde_, 0, 0, RSQrqt_tilde_, 0, 0);
    }
    diaad(nx, 1.0, D_x, 0, RSQrqt_tilde_, 0, 0);
    // For the regularized branch there is no eq elimination, so Ppt_ tracks RSQrqt_tilde_.
    gecp(nx + 1, nx, RSQrqt_tilde_, 0, 0, Ppt_, 0, 0);

    potrf_l_mn(nx + 1, nx, Ppt_, 0, 0, Llt_, 0, 0);
    if (!check_reg(nx, &Llt_.mat(), 0, 0))
        return LinsolReturnFlag::INDEFINITE;
    rank_eq_ = 0;

    // ----- 2. Forward substitute the primal direction.
    rowex(nx, -1.0, Llt_, nx, 0, x, 0);
    trsv_ltn(nx, Llt_, 0, 0, x, 0, x, 0);

    // ----- 3. Recover the multipliers from the (now eliminated) penalty rows.
    if (ng > 0)
    {
        gemv_t(nx, ng, 1.0, jacobian.Gg_eqt, 0, 0, x, 0, 1.0, g,
               0, eq_mult, 0);
        eq_mult.block(ng, 0) =
            eq_mult.block(ng, 0) / D_eq.block(ng, 0);
    }
    if (ng_ineq > 0)
    {
        gemv_t(nx, ng_ineq, 1.0, jacobian.Gg_ineqt, 0, 0, x, 0, 1.0, g,
               info.offset_g_eq_slack, eq_mult, info.offset_g_eq_slack);
        eq_mult.block(ng_ineq, info.offset_g_eq_slack) =
            eq_mult.block(ng_ineq, info.offset_g_eq_slack) / D_s.block(ng_ineq, 0);
    }
    return LinsolReturnFlag::SUCCESS;
}

LinsolReturnFlag AugSystemSolver<DenseType>::solve_rhs(const ProblemInfo<DenseType> &info,
                                                       const Jacobian<DenseType> &jacobian,
                                                       const Hessian<DenseType> &hessian,
                                                       const VecRealView &D_s,
                                                       const VecRealView &f, const VecRealView &g,
                                                       VecRealView &x, VecRealView &eq_mult)
{
    const Index nx = info.dims.nx_tangent;
    const Index ng = info.dims.ng;
    const Index ng_ineq = info.dims.ng_ineq;

    // Build the right-hand side of the reduced system (only RHS vector ops, no refactor).
    veccp(nx, f, 0, v_RSQrqt_tilde_, 0);
    if (ng_ineq > 0)
    {
        for (Index i = 0; i < ng_ineq; i++)
        {
            const Scalar di = D_s(i);
            const Scalar gi = g(info.offset_g_eq_slack + i);
            v_Ggt_ineq_(i) = gi / di;
        }
        gemv_n(nx, ng_ineq, 1.0, jacobian.Gg_ineqt, 0, 0, v_Ggt_ineq_, 0, 1.0, v_RSQrqt_tilde_, 0,
               v_RSQrqt_tilde_, 0);
    }

    if (ng > 0)
    {
        // v_Ggt_stripe = Pl @ g_eq, then apply U⁻¹ from the LU factors stored in Ggt_stripe_.
        veccp(ng, g, 0, v_Ggt_stripe_, 0);
        Pl_.apply(rank_eq_, &v_Ggt_stripe_.vec(), 0);
        trsv_utu(rank_eq_, Ggt_stripe_, 0, 0, v_Ggt_stripe_, 0, v_Ggt_stripe_, 0);

        // v_Ggt_tilde = -L⁻ᵀ v_Ggt_stripe.
        veccpsc(rank_eq_, -1.0, v_Ggt_stripe_, 0, v_Ggt_tilde_, 0);
        trsv_ltn(rank_eq_, Ggt_stripe_, 0, 0, v_Ggt_tilde_, 0, v_Ggt_tilde_, 0);

        // Permute v_RSQrqt_tilde to the (Pr-permuted) Hessian's coordinates. After this
        // line v_RSQrqt_tilde_ plays the role of the "v_Ppt" vector — we keep it
        // untouched until the equality-multiplier recovery below.
        Pr_.apply(rank_eq_, &v_RSQrqt_tilde_.vec(), 0);

        // v_GgLt = v_RSQrqt_tilde + Ppt[:nx, :rank] @ v_Ggt_tilde
        veccp(nx, v_RSQrqt_tilde_, 0, v_GgLt_, 0);
        gemv_n(nx, rank_eq_, 1.0, Ppt_, 0, 0, v_Ggt_tilde_, 0, 1.0, v_GgLt_, 0, v_GgLt_, 0);
        // Schur-eliminate: v_RSQrqt_hat = v_GgLt[rank:] + Ggt_tilde @ v_GgLt[:rank]
        gemv_n(nx - rank_eq_, rank_eq_, 1.0, Ggt_tilde_, 0, 0, v_GgLt_, 0, 1.0, v_GgLt_, rank_eq_,
               v_RSQrqt_hat_, 0);
        trsv_lnn(nx - rank_eq_, Llt_, 0, 0, v_RSQrqt_hat_, 0, v_Llt_, 0);

        // Recover x and eq_mult by forward substitution (mirror of solve()).
        veccpsc(nx - rank_eq_, -1.0, v_Llt_, 0, x, rank_eq_);
        trsv_ltn(nx - rank_eq_, Llt_, 0, 0, x, rank_eq_, x,
                 rank_eq_);
        veccp(rank_eq_, v_Ggt_tilde_, 0, x, 0);
        gemv_t(nx - rank_eq_, rank_eq_, 1.0, Ggt_tilde_, 0, 0, x,
               rank_eq_, 1.0, x, 0, x,
               0);
        veccpsc(rank_eq_, -1.0, v_RSQrqt_tilde_, 0, eq_mult, 0);
        gemv_t(nx, rank_eq_, -1.0, Ppt_, 0, 0, x, 0, 1.0, eq_mult,
               0, eq_mult, 0);
        trsv_lnn(rank_eq_, Ggt_stripe_, 0, 0, eq_mult, 0, eq_mult,
                 0);
        trsv_unu(rank_eq_, rank_eq_, Ggt_stripe_, 0, 0, eq_mult, 0, eq_mult,
                 0);
        Pl_.apply_inverse(rank_eq_, &eq_mult.vec(), 0);
        Pr_.apply_inverse(rank_eq_, &x.vec(), 0);
    }
    else
    {
        trsv_lnn(nx, Llt_, 0, 0, v_RSQrqt_tilde_, 0, v_Llt_, 0);
        veccpsc(nx, -1.0, v_Llt_, 0, x, 0);
        trsv_ltn(nx, Llt_, 0, 0, x, 0, x, 0);
    }

    if (ng_ineq > 0)
    {
        gemv_t(nx, ng_ineq, 1.0, jacobian.Gg_ineqt, 0, 0, x, 0, 1.0, g,
               info.offset_g_eq_slack, eq_mult, info.offset_g_eq_slack);
        eq_mult.block(ng_ineq, info.offset_g_eq_slack) =
            eq_mult.block(ng_ineq, info.offset_g_eq_slack) / D_s.block(ng_ineq, 0);
    }
    return LinsolReturnFlag::SUCCESS;
}

LinsolReturnFlag AugSystemSolver<DenseType>::solve_rhs(const ProblemInfo<DenseType> &info,
                                                       const Jacobian<DenseType> &jacobian,
                                                       const Hessian<DenseType> &hessian,
                                                       const VecRealView &D_eq,
                                                       const VecRealView &D_s,
                                                       const VecRealView &f, const VecRealView &g,
                                                       VecRealView &x, VecRealView &eq_mult)
{
    const Index nx = info.dims.nx_tangent;
    const Index ng = info.dims.ng;
    const Index ng_ineq = info.dims.ng_ineq;

    // Form the reduced RHS, reusing the Cholesky factor produced by the regularized solve().
    veccp(nx, f, 0, v_RSQrqt_tilde_, 0);
    if (ng > 0)
    {
        for (Index i = 0; i < ng; i++)
            v_Ggt_stripe_(i) = g(0 + i) / D_eq(i);
        gemv_n(nx, ng, 1.0, jacobian.Gg_eqt, 0, 0, v_Ggt_stripe_, 0, 1.0, v_RSQrqt_tilde_, 0,
               v_RSQrqt_tilde_, 0);
    }
    if (ng_ineq > 0)
    {
        for (Index i = 0; i < ng_ineq; i++)
            v_Ggt_ineq_(i) = g(info.offset_g_eq_slack + i) / D_s(i);
        gemv_n(nx, ng_ineq, 1.0, jacobian.Gg_ineqt, 0, 0, v_Ggt_ineq_, 0, 1.0, v_RSQrqt_tilde_, 0,
               v_RSQrqt_tilde_, 0);
    }
    trsv_lnn(nx, Llt_, 0, 0, v_RSQrqt_tilde_, 0, v_Llt_, 0);
    veccpsc(nx, -1.0, v_Llt_, 0, x, 0);
    trsv_ltn(nx, Llt_, 0, 0, x, 0, x, 0);

    if (ng > 0)
    {
        gemv_t(nx, ng, 1.0, jacobian.Gg_eqt, 0, 0, x, 0, 1.0, g,
               0, eq_mult, 0);
        eq_mult.block(ng, 0) =
            eq_mult.block(ng, 0) / D_eq.block(ng, 0);
    }
    if (ng_ineq > 0)
    {
        gemv_t(nx, ng_ineq, 1.0, jacobian.Gg_ineqt, 0, 0, x, 0, 1.0, g,
               info.offset_g_eq_slack, eq_mult, info.offset_g_eq_slack);
        eq_mult.block(ng_ineq, info.offset_g_eq_slack) =
            eq_mult.block(ng_ineq, info.offset_g_eq_slack) / D_s.block(ng_ineq, 0);
    }
    return LinsolReturnFlag::SUCCESS;
}
