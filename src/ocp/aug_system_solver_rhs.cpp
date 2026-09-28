// Copyright (c) 2026 Lander Vanroye, KU Leuven
// Matrix-RHS overloads reuse the native augmented-system factorization.
#include "fatrop/ocp/aug_system_solver.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include "fatrop/linear_algebra/batched_rhs.hpp"
#include <optional>

using namespace fatrop;

void AugSystemSolver<OcpType>::allocate_rhs_workspace()
{
    batched::allocate_like(rhs_many_work_[0], v_Ppt, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[1], v_Hh, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[2], v_AL, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[3], v_RSQrqt_tilde, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[4], v_Ggt_stripe, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[5], v_Ggt_tilde, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[6], v_GgLt, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[7], v_RSQrqt_hat, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[8], v_Llt, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[9], v_Llt_shift, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[10], v_GgIt_tilde, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[11], v_GgLIt, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[12], v_HhIt, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[13], v_PpIt_hat, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[14], v_LlIt, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[15], v_Ggt_ineq_temp, rhs_capacity_);
    batched::allocate_like(rhs_many_work_[16], v_tmp, rhs_capacity_);
}

LinsolReturnFlag AugSystemSolver<OcpType>::solve_rhs(const ProblemInfo<OcpType> &info,
                                               const Jacobian<OcpType> &jacobian,
                                               const Hessian<OcpType> &hessian,
                                               const VecRealView &D_s, const MatRealView &f,
                                               const MatRealView &g, MatRealView &x,
                                               MatRealView &eq_mult)
{
    batched::check_dimensions(info.number_of_tangent_variables, info.number_of_eq_constraints,
                              f, g, x, eq_mult);
    if (f.n() == 0) return LinsolReturnFlag::SUCCESS;
    if (f.n() > rhs_capacity_)
        return batched::solve_panels(rhs_capacity_, f, g, x, eq_mult,
            [&](const MatRealView &fp, const MatRealView &gp, MatRealView &xp, MatRealView &lp)
            { return solve_rhs(info, jacobian, hessian, D_s, fp, gp, xp, lp); });
    batched::WorkspaceView m_Ppt(rhs_many_work_[0], f.n());
    batched::WorkspaceView m_Hh(rhs_many_work_[1], f.n());
    batched::WorkspaceView m_AL(rhs_many_work_[2], f.n());
    batched::WorkspaceView m_RSQrqt_tilde(rhs_many_work_[3], f.n());
    batched::WorkspaceView m_Ggt_stripe(rhs_many_work_[4], f.n());
    batched::WorkspaceView m_Ggt_tilde(rhs_many_work_[5], f.n());
    batched::WorkspaceView m_GgLt(rhs_many_work_[6], f.n());
    batched::WorkspaceView m_RSQrqt_hat(rhs_many_work_[7], f.n());
    batched::WorkspaceView m_Llt(rhs_many_work_[8], f.n());
    batched::WorkspaceView m_Llt_shift(rhs_many_work_[9], f.n());
    batched::WorkspaceView m_GgIt_tilde(rhs_many_work_[10], f.n());
    batched::WorkspaceView m_GgLIt(rhs_many_work_[11], f.n());
    batched::WorkspaceView m_HhIt(rhs_many_work_[12], f.n());
    batched::WorkspaceView m_PpIt_hat(rhs_many_work_[13], f.n());
    batched::WorkspaceView m_LlIt(rhs_many_work_[14], f.n());
    batched::WorkspaceView m_Ggt_ineq_temp(rhs_many_work_[15], f.n());
    batched::WorkspaceView m_tmp(rhs_many_work_[16], f.n());

    std::optional<MatRealView> m_RSQrq_hat_curr_p;
    Index rank_k;

    // Backward recursion: eliminate controls and propagate constraint RHS blocks.
    for (Index k = info.dims.K - 1; k >= 0; --k)
    {
        const Index nu = info.dims.number_of_tangent_controls[k];
        const Index nx = info.dims.number_of_tangent_states[k];
        const Index ng = info.dims.number_of_eq_constraints[k];
        const Index ng_ineq = info.dims.number_of_ineq_constraints[k];
        const Index offset_ineq_k = info.offsets_slack[k];
        const Index offs_g_ineq_k = info.offsets_g_eq_slack[k];
        const Index offs_g_k = info.offsets_g_eq_path[k];
        const Index offs = info.offsets_tangent_u[k];

        Index gamma_k;
        if (k == info.dims.K - 1)
        {
            gamma_k = ng;
            gamma[k] = gamma_k;
            batched::gecp(ng, g, offs_g_k, m_Ggt_stripe[0], 0);
            batched::gecp(nu + nx, f, offs, m_RSQrqt_tilde[k], 0);
        }
        else
        {
            const Index offs_dyn_k = info.offsets_g_eq_dyn[k];
            const Index nxp1 = info.dims.number_of_tangent_states[k + 1];
            const Index Hp1_size = gamma[k + 1] - rho[k + 1];
            gamma_k = Hp1_size + ng;
            batched::gemm_nn(nxp1, nxp1, 1.0, Ppt[k + 1], 0, 0, g, offs_dyn_k, 0.0, m_AL[0], 0, m_AL[0], 0);
            batched::gead(nxp1, 1.0, m_Ppt[k + 1], 0, m_AL[0], 0, m_AL[0], 0);
            batched::gemm_nn(nu + nx, nxp1, 1.0, jacobian.BAbt[k], 0, 0, m_AL[0], 0, 1.0, f, offs,
                   m_RSQrqt_tilde[k], 0);
            if (gamma_k > 0)
            {
                if (ng > 0)
                {
                    batched::gecp(ng, g, offs_g_k, m_Ggt_stripe[0], 0);
                }
                if (Hp1_size > 0)
                {
                    batched::gemm_nn(Hp1_size, nxp1, 1.0, Hh[k + 1], 0, 0, g, offs_dyn_k, 0.0,
                           m_Ggt_stripe[0], ng, m_Ggt_stripe[0], ng);
                    batched::gead(Hp1_size, 1.0, m_Hh[k + 1], 0, m_Ggt_stripe[0], ng, m_Ggt_stripe[0], ng);
                }
            }
            else
            {
                rank_k = 0;
                m_RSQrq_hat_curr_p.emplace(m_RSQrqt_tilde[k]);
            }
        }
        if (ng_ineq > 0)
        {
            batched::divide_rows(ng_ineq, g, offs_g_ineq_k, D_s, offset_ineq_k, m_Ggt_ineq_temp[0], 0);
            batched::gemm_nn(nu + nx, ng_ineq, 1.0, jacobian.Gg_ineqt[k], 0, 0, m_Ggt_ineq_temp[0], 0, 1.0,
                   m_RSQrqt_tilde[k], 0, m_RSQrqt_tilde[k], 0);
        }
        {
            rank_k = rho[k];
            gecp(rank_k, gamma_k, Ggt_tilde[k], nu - rank_k + nx + 1, 0, Ggt_stripe[0], 0, 0);
            batched::rowpe(Pl[k], rank_k, m_Ggt_stripe[0], 0);
            batched::trsm_lutu(rank_k, Ggt_stripe[0], 0, 0, m_Ggt_stripe[0], 0, m_Ggt_stripe[0], 0);
            batched::gemm_tn(rank_k, gamma_k - rank_k, -1.0, Ggt_stripe[0], 0, rank_k, m_Ggt_stripe[0], 0,
                   1.0, m_Ggt_stripe[0], rank_k, m_Ggt_stripe[0], rank_k);

            if (gamma_k - rank_k > 0)
            {
                batched::gecp(gamma_k - rank_k, m_Ggt_stripe[0], rank_k, m_Hh[k], 0);
            }
            if (rank_k > 0)
            {
                batched::gecpsc(rank_k, -1.0, m_Ggt_stripe[0], 0, m_Ggt_tilde[k], 0);
                batched::trsm_lltn(rank_k, Ggt_stripe[0], 0, 0, m_Ggt_tilde[k], 0, m_Ggt_tilde[k], 0);
                batched::rowpe(Pr[k], rank_k, m_RSQrqt_tilde[k], 0);
                batched::gecp(nu + nx, m_RSQrqt_tilde[k], 0, m_GgLt[0], 0);
                batched::gemm_nn(nu + nx, rank_k, 1.0, RSQrqt_tilde[k], 0, 0, m_Ggt_tilde[k], 0, 1.0,
                       m_GgLt[0], 0, m_GgLt[0], 0);
                batched::gemm_nn(nu + nx - rank_k, rank_k, 1.0, Ggt_tilde[k], 0, 0, m_GgLt[0], 0, 1.0,
                       m_GgLt[0], rank_k, m_RSQrqt_hat[0], 0);
                m_RSQrq_hat_curr_p.emplace(m_RSQrqt_hat[0]);
            }
            else
            {
                m_RSQrq_hat_curr_p.emplace(m_RSQrqt_tilde[k]);
            }
        }

        {
            if (nu - rank_k > 0)
            {
                batched::trsm_llnn(nu - rank_k, Llt[k], 0, 0, *m_RSQrq_hat_curr_p, 0, m_Llt[k], 0);
                gecp(nx + 1, nu - rank_k, Llt[k], nu - rank_k, 0, Llt_shift[0], 0, 0);
                batched::gecp(nu - rank_k, m_Llt[k], 0, m_Llt_shift[0], 0);
                batched::gecp(nx, *m_RSQrq_hat_curr_p, nu - rank_k, m_Ppt[k], 0);
                batched::gemm_nn(nx, nu - rank_k, -1.0, Llt_shift[0], 0, 0, m_Llt_shift[0], 0, 1.0, m_Ppt[k],
                       0, m_Ppt[k], 0);
                if (increased_accuracy)
                {
                    batched::gemm_tn(nu - rank_k, gamma_k - rank_k, -1.0, Ggt_tilde[k], 0, rank_k, m_Llt[k],
                           0, 1.0, m_Hh[k], 0, m_Hh[k], 0);
                }
            }
            else
            {
                batched::gecp(nx, *m_RSQrq_hat_curr_p, 0, m_Ppt[k], 0);
            }
        }
    }
    {
        const Index nx = info.dims.number_of_tangent_states[0];
        Index gamma_I = gamma[0] - rho[0];
        if (gamma_I > 0)
        {
            batched::gecp(gamma_I, m_Hh[0], 0, m_HhIt[0], 0);
            batched::rowpe(PlI[0], rankI, m_HhIt[0], 0);
            batched::trsm_lutu(rankI, HhIt[0], 0, 0, m_HhIt[0], 0, m_HhIt[0], 0);
            batched::gemm_tn(rankI, gamma_I - rankI, -1.0, HhIt[0], 0, rankI, m_HhIt[0], 0, 1.0, m_HhIt[0],
                   rankI, m_HhIt[0], rankI);
            batched::gecpsc(rankI, -1.0, m_HhIt[0], 0, m_GgIt_tilde[0], 0);
            batched::trsm_lltn(rankI, HhIt[0], 0, 0, m_GgIt_tilde[0], 0, m_GgIt_tilde[0], 0);
            batched::rowpe(PrI[0], rankI, m_Ppt[0], 0);
            batched::gecp(nx, m_Ppt[0], 0, m_GgLIt[0], 0);
            batched::gemm_nn(nx, rankI, 1.0, Ppt[0], 0, 0, m_GgIt_tilde[0], 0, 1.0, m_GgLIt[0], 0, m_GgLIt[0],
                   0);
            batched::gemm_nn(nx - rankI, rankI, 1.0, GgIt_tilde[0], 0, 0, m_GgLIt[0], 0, 1.0, m_GgLIt[0],
                   rankI, m_PpIt_hat[0], 0);
            batched::trsm_llnn(nx - rankI, LlIt[0], 0, 0, m_PpIt_hat[0], 0, m_LlIt[0], 0);
        }
        else
        {
            batched::trsm_llnn(nx, LlIt[0], 0, 0, m_Ppt[0], 0, m_LlIt[0], 0);
        }
    }
    {
        const Index nx = info.dims.number_of_tangent_states[0];
        const Index nu = info.dims.number_of_tangent_controls[0];
        const Index offs_u = info.offsets_tangent_u[0];
        const Index offs_x = info.offsets_tangent_x[0];
        const Index offs_g = info.offsets_g_eq_path[0];
        batched::gecpsc(nx - rankI, -1.0, m_LlIt[0], 0, x, offs_x + rankI);
        batched::trsm_lltn(nx - rankI, LlIt[0], 0, 0, x, offs_x + rankI, x, offs_x + rankI);
        batched::gecp(rankI, m_GgIt_tilde[0], 0, x, offs_x);
        batched::gemm_tn(nx - rankI, rankI, 1.0, GgIt_tilde[0], 0, 0, x, offs_x + rankI, 1.0, x, offs_x, x,
               offs_x);
        batched::gecpsc(rankI, -1.0, m_Ppt[0], 0, eq_mult, offs_g);
        batched::gemm_tn(nx, rankI, -1.0, Ppt[0], 0, 0, x, nu, 1.0, eq_mult, offs_g, eq_mult, offs_g);
        batched::trsm_llnn(rankI, HhIt[0], 0, 0, eq_mult, offs_g, eq_mult, offs_g);
        batched::trsm_lunu(rankI, rankI, HhIt[0], 0, 0, eq_mult, offs_g, eq_mult, offs_g);
        batched::rowpei(PlI[0], rankI, eq_mult, offs_g);
        batched::rowpei(PrI[0], rankI, x, offs_x);
    }
    // Forward recursion and multiplier recovery.
    for (Index k = 0; k < info.dims.K; k++)
    {

        const Index nx = info.dims.number_of_tangent_states[k];
        const Index nu = info.dims.number_of_tangent_controls[k];
        const Index offs = info.offsets_tangent_u[k];
        const Index offs_x = info.offsets_tangent_x[k];
        const Index rho_k = rho[k];
        const Index numrho_k = nu - rho_k;
        const Index offs_g_k = info.offsets_g_eq_path[k];
        const Index gammamrho_k = gamma[k] - rho[k];
        const Index gamma_k = gamma[k];
        const Index ng_ineq = info.dims.number_of_ineq_constraints[k];
        const Index offs_eq_ineq = info.offsets_g_eq_slack[k];
        const Index offs_slack = info.offsets_slack[k];
        if (numrho_k > 0)
        {
            batched::gecpsc(numrho_k, -1.0, m_Llt[k], 0, x, offs + rho_k);
            if (increased_accuracy)
            {
                batched::gemm_nn(nu - rho_k, gamma_k - rho_k, -1.0, Ggt_tilde[k], 0, rho_k, eq_mult, offs_g_k,
                       1.0, x, offs + rho_k, x, offs + rho_k);
            }
            batched::gemm_tn(nx, numrho_k, -1.0, Llt[k], numrho_k, 0, x, offs_x, 1.0, x, offs + rho_k, x,
                   offs + rho_k);
            batched::trsm_lltn(numrho_k, Llt[k], 0, 0, x, offs + rho_k, x, offs + rho_k);
        }

        if (rho_k > 0)
        {

            batched::gecp(rho_k, m_Ggt_tilde[k], 0, x, offs);
            batched::gemm_tn(nx + numrho_k, rho_k, 1.0, Ggt_tilde[k], 0, 0, x, offs + rho_k, 1.0, x, offs, x,
                   offs);
            batched::gecp(gammamrho_k, eq_mult, offs_g_k, m_tmp[0], 0);
            batched::gecp(gammamrho_k, m_tmp[0], 0, eq_mult, offs_g_k + rho_k);
            batched::gecpsc(rho_k, -1.0, m_RSQrqt_tilde[k], 0, eq_mult, offs_g_k);
            batched::gemm_tn(nu + nx, rho_k, -1.0, RSQrqt_tilde[k], 0, 0, x, offs, 1.0, eq_mult, offs_g_k,
                   eq_mult, offs_g_k);
            gecp(rho_k, gamma_k, Ggt_tilde[k], nu - rho_k + nx + 1, 0, AL[0], 0, 0);
            batched::trsm_llnn(rho_k, AL[0], 0, 0, eq_mult, offs_g_k, eq_mult, offs_g_k);
            batched::trsm_lunu(rho_k, gamma_k, AL[0], 0, 0, eq_mult, offs_g_k, eq_mult, offs_g_k);
            batched::rowpei(Pl[k], rho_k, eq_mult, offs_g_k);
            batched::rowpei(Pr[k], rho_k, x, offs);
        }
        if (ng_ineq > 0)
        {
            batched::gemm_tn(nu + nx, ng_ineq, 1.0, jacobian.Gg_ineqt[k], 0, 0, x, offs, 1.0, g, offs_eq_ineq,
                   eq_mult, offs_eq_ineq);
            batched::divide_rows(ng_ineq, eq_mult, offs_eq_ineq, D_s, offs_slack, eq_mult, offs_eq_ineq);
        }
        if (k != info.dims.K - 1)
        {
            const Index nxp1 = info.dims.number_of_tangent_states[k + 1];
            const Index nup1 = info.dims.number_of_tangent_controls[k + 1];
            const Index offsp1 = info.offsets_tangent_u[k + 1];
            const Index offsxp1 = info.offsets_tangent_x[k + 1];
            const Index offs_g_kp1 = info.offsets_g_eq_path[k + 1];
            const Index offs_dyn_k = info.offsets_g_eq_dyn[k];
            const Index gammamrho_kp1 = gamma[k + 1] - rho[k + 1];
            batched::gecp(nxp1, g, offs_dyn_k, x, offsxp1);
            batched::gemm_tn(nu + nx, nxp1, 1.0, jacobian.BAbt[k], 0, 0, x, offs, 1.0, x, offsxp1, x,
                   offsxp1);
            batched::gecp(nxp1, m_Ppt[k + 1], 0, eq_mult, offs_dyn_k);
            batched::gemm_tn(nxp1, nxp1, 1.0, Ppt[k + 1], 0, 0, x, offsxp1, 1.0, eq_mult, offs_dyn_k, eq_mult,
                   offs_dyn_k);
            batched::gemm_tn(gammamrho_kp1, nxp1, 1.0, Hh[k + 1], 0, 0, eq_mult, offs_g_kp1, 1.0, eq_mult,
                   offs_dyn_k, eq_mult, offs_dyn_k);
        }
    }
    return LinsolReturnFlag::SUCCESS;
}
LinsolReturnFlag AugSystemSolver<OcpType>::solve_rhs(const ProblemInfo<OcpType> &info,
                                               const Jacobian<OcpType> &jacobian,
                                               const Hessian<OcpType> &hessian,
                                               const VecRealView &D_eq, const VecRealView &D_s,
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
    batched::WorkspaceView m_Ppt(rhs_many_work_[0], f.n());
    batched::WorkspaceView m_Hh(rhs_many_work_[1], f.n());
    batched::WorkspaceView m_AL(rhs_many_work_[2], f.n());
    batched::WorkspaceView m_RSQrqt_tilde(rhs_many_work_[3], f.n());
    batched::WorkspaceView m_Ggt_stripe(rhs_many_work_[4], f.n());
    batched::WorkspaceView m_Ggt_tilde(rhs_many_work_[5], f.n());
    batched::WorkspaceView m_GgLt(rhs_many_work_[6], f.n());
    batched::WorkspaceView m_RSQrqt_hat(rhs_many_work_[7], f.n());
    batched::WorkspaceView m_Llt(rhs_many_work_[8], f.n());
    batched::WorkspaceView m_Llt_shift(rhs_many_work_[9], f.n());
    batched::WorkspaceView m_GgIt_tilde(rhs_many_work_[10], f.n());
    batched::WorkspaceView m_GgLIt(rhs_many_work_[11], f.n());
    batched::WorkspaceView m_HhIt(rhs_many_work_[12], f.n());
    batched::WorkspaceView m_PpIt_hat(rhs_many_work_[13], f.n());
    batched::WorkspaceView m_LlIt(rhs_many_work_[14], f.n());
    batched::WorkspaceView m_Ggt_ineq_temp(rhs_many_work_[15], f.n());
    batched::WorkspaceView m_tmp(rhs_many_work_[16], f.n());

    std::optional<MatRealView> m_RSQrq_hat_curr_p;
    // Backward recursion: eliminate controls and propagate constraint RHS blocks.
    for (Index k = info.dims.K - 1; k >= 0; --k)
    {
        const Index offs_ux_k = info.offsets_tangent_u[k];
        const Index nu = info.dims.number_of_tangent_controls[k];
        const Index nx = info.dims.number_of_tangent_states[k];
        const Index ng = info.dims.number_of_eq_constraints[k];
        const Index ng_ineq = info.dims.number_of_ineq_constraints[k];
        const Index offs_g_dyn = info.offsets_g_eq_dyn[k];
        const Index offs_g_eq = info.offsets_g_eq_path[k];
        const Index offs_ge_eq_ineq = info.offsets_g_eq_slack[k];

        if (k == info.dims.K - 1)
        {
            batched::gecp(nu + nx, f, offs_ux_k, m_RSQrqt_tilde[k], 0);
        }
        else
        {
            const Index nxp1 = info.dims.number_of_tangent_states[k + 1];
            batched::gemm_nn(nxp1, nxp1, 1.0, Ppt[k + 1], 0, 0, g, offs_g_dyn, 0.0, m_AL[0], 0, m_AL[0], 0);
            batched::gead(nxp1, 1.0, m_Ppt[k + 1], 0, m_AL[0], 0, m_AL[0], 0);
            batched::gemm_nn(nu + nx, nxp1, 1.0, jacobian.BAbt[k], 0, 0, m_AL[0], 0, 1.0, f, offs_ux_k,
                   m_RSQrqt_tilde[k], 0);
        }
        if (ng > 0)
        {
            const Index offs_eq_k = info.offsets_eq[k];
            batched::divide_rows(ng, g, offs_g_eq, D_eq, offs_eq_k, m_Ggt_stripe[0], 0);
            batched::gemm_nn(nu + nx, ng, 1.0, jacobian.Gg_eqt[k], 0, 0, m_Ggt_stripe[0], 0, 1.0,
                   m_RSQrqt_tilde[k], 0, m_RSQrqt_tilde[k], 0);
        }
        if (ng_ineq > 0)
        {
            const Index offs_ineq_k = info.offsets_slack[k];
            batched::divide_rows(ng_ineq, g, offs_ge_eq_ineq, D_s, offs_ineq_k, m_Ggt_ineq_temp[0], 0);
            batched::gemm_nn(nu + nx, ng_ineq, 1.0, jacobian.Gg_ineqt[k], 0, 0, m_Ggt_ineq_temp[0], 0, 1.0,
                   m_RSQrqt_tilde[k], 0, m_RSQrqt_tilde[k], 0);
        }
        {
            m_RSQrq_hat_curr_p.emplace(m_RSQrqt_tilde[k]);
        }
        {
            batched::trsm_llnn(nu, Llt[k], 0, 0, *m_RSQrq_hat_curr_p, 0, m_Llt[k], 0);
            batched::gecp(nu, m_Llt[k], 0, m_Llt_shift[0], 0);
            batched::gemm_nn(nx, nu, -1.0, Llt[k], nu, 0, m_Llt_shift[0], 0, 1.0, m_RSQrqt_tilde[k], nu,
                   m_Ppt[k], 0);
        }
    }
    {
        const Index nx = info.dims.number_of_tangent_states[0];
        {
            batched::trsm_llnn(nx, LlIt[0], 0, 0, m_Ppt[0], 0, m_LlIt[0], 0);
        }
    }
    {
        const Index nx = info.dims.number_of_tangent_states[0];
        const Index nu = info.dims.number_of_tangent_controls[0];
        const Index offs_x = info.offsets_tangent_x[0];
        batched::gecpsc(nx, -1.0, m_LlIt[0], 0, x, offs_x);
        batched::trsm_lltn(nx, LlIt[0], 0, 0, x, offs_x, x, offs_x);
    }
    for (Index k = 0; k < info.dims.K; k++)
    {
        const Index nx = info.dims.number_of_tangent_states[k];
        const Index nu = info.dims.number_of_tangent_controls[k];
        const Index offs = info.offsets_tangent_u[k];
        const Index offs_x = info.offsets_tangent_x[k];
        const Index offs_dyn_eq_k = info.offsets_g_eq_dyn[k];
        batched::gecpsc(nu, -1.0, m_Llt[k], 0, x, offs);
        batched::gemm_tn(nx, nu, -1.0, Llt[k], nu, 0, x, offs_x, 1.0, x, offs, x, offs);
        batched::trsm_lltn(nu, Llt[k], 0, 0, x, offs, x, offs);
        if (k != info.dims.K - 1)
        {
            const Index nxp1 = info.dims.number_of_tangent_states[k + 1];
            const Index offsp1 = info.offsets_tangent_u[k + 1];
            const Index offs_x_p1 = info.offsets_tangent_x[k + 1];
            batched::gecp(nxp1, g, offs_dyn_eq_k, x, offs_x_p1);
            batched::gemm_tn(nu + nx, nxp1, 1.0, jacobian.BAbt[k], 0, 0, x, offs, 1.0, x, offs_x_p1, x,
                   offs_x_p1);
            batched::gecp(nxp1, m_Ppt[k + 1], 0, eq_mult, offs_dyn_eq_k);
            batched::gemm_tn(nxp1, nxp1, 1.0, Ppt[k + 1], 0, 0, x, offs_x_p1, 1.0, eq_mult,
                   offs_dyn_eq_k, eq_mult, offs_dyn_eq_k);
        }
    }

    for (Index k = 0; k < info.dims.K; k++)
    {
        const Index nx = info.dims.number_of_tangent_states[k];
        const Index nu = info.dims.number_of_tangent_controls[k];
        const Index ng = info.dims.number_of_eq_constraints[k];
        const Index offs = info.offsets_tangent_u[k];
        const Index offs_g_k = info.offsets_g_eq_path[k];
        const Index offs_eq = info.offsets_eq[k];
        if (ng > 0)
        {
            batched::gemm_tn(nu + nx, ng, 1.0, jacobian.Gg_eqt[k], 0, 0, x, offs, 1.0, g, offs_g_k,
                   eq_mult, offs_g_k);
            batched::divide_rows(ng, eq_mult, offs_g_k, D_eq, offs_eq, eq_mult, offs_g_k);
        }
    }

    for (Index k = 0; k < info.dims.K; k++)
    {
        const Index nx = info.dims.number_of_tangent_states[k];
        const Index nu = info.dims.number_of_tangent_controls[k];
        const Index ng_ineq = info.dims.number_of_ineq_constraints[k];
        const Index offs = info.offsets_tangent_u[k];
        const Index offs_gineq_k = info.offsets_g_eq_slack[k];
        const Index offs_slack = info.offsets_slack[k];
        if (ng_ineq > 0)
        {
            batched::gemm_tn(nu + nx, ng_ineq, 1.0, jacobian.Gg_ineqt[k], 0, 0, x, offs, 1.0, g, offs_gineq_k,
                   eq_mult, offs_gineq_k);
            batched::divide_rows(ng_ineq, eq_mult, offs_gineq_k, D_s, offs_slack, eq_mult, offs_gineq_k);
        }
    }
    return LinsolReturnFlag::SUCCESS;
}