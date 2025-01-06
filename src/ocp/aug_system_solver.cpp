//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#include "fatrop/ocp/aug_system_solver.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include <algorithm>
using namespace fatrop;

bool check_reg(const Index m, MAT *sA, const Index ai, const Index aj)
{
    for (Index i = 0; i < m; i++)
    {
        if (blasfeo_matel_wrap(sA, ai + i, aj + i) < 1e-8)
            return false;
    }
    return true;
}

OcpAugSystemSolver::OcpAugSystemSolver(const ProblemInfo<OcpType> &info)
{
    Index max_number_of_controls =
        *std::max_element(info.dims.number_of_controls.begin(), info.dims.number_of_controls.end());
    Index max_number_of_states =
        *std::max_element(info.dims.number_of_states.begin(), info.dims.number_of_states.end());
    Index max_number_of_variables = *std::max_element(info.number_of_stage_variables.begin(),
                                                      info.number_of_stage_variables.end());
    Index max_number_of_ineq_constraints = *std::max_element(
        info.dims.number_of_ineq_constraints.begin(), info.dims.number_of_ineq_constraints.end());
    Index max_number_of_eq_consttraints = *std::max_element(
        info.dims.number_of_eq_constraints.begin(), info.dims.number_of_eq_constraints.end());

    AL.emplace_back(max_number_of_variables + 1, max_number_of_variables);
    Ggt_stripe.emplace_back(max_number_of_variables + 1, max_number_of_variables);
    GgLt.emplace_back(max_number_of_variables + 1, max_number_of_variables);
    RSQrqt_hat.emplace_back(max_number_of_variables + 1, max_number_of_variables);
    Llt_shift.emplace_back(max_number_of_variables + 1, max_number_of_controls);
    GgIt_tilde.emplace_back(info.dims.number_of_states[0] + 1, info.dims.number_of_states[0]);
    GgLIt.emplace_back(info.dims.number_of_states[0] + 1, info.dims.number_of_states[0]);
    HhIt.emplace_back(info.dims.number_of_states[0] + 1, info.dims.number_of_states[0]);
    PpIt_hat.emplace_back(info.dims.number_of_states[0] + 1, info.dims.number_of_states[0]);
    LlIt.emplace_back(info.dims.number_of_states[0] + 1, info.dims.number_of_states[0]);
    Ggt_ineq_temp.emplace_back(max_number_of_variables + 1, max_number_of_ineq_constraints);

    Ppt.reserve(info.dims.K);
    Hh.reserve(info.dims.K);
    RSQrqt_tilde.reserve(info.dims.K);
    Ggt_tilde.reserve(info.dims.K);
    Llt.reserve(info.dims.K);
    for (Index k = 0; k < info.dims.K; k++)
    {
        Index nu = info.dims.number_of_controls[k];
        Index nx = info.dims.number_of_states[k];
        Index ng_ineq = info.dims.number_of_ineq_constraints[k];
        Index ng_eq = info.dims.number_of_eq_constraints[k];
        Ppt.emplace_back(nx + 1, nx);
        Hh.emplace_back(nx, nx + 1);
        RSQrqt_tilde.emplace_back(nu + nx + 1, nx + nu);
        Ggt_tilde.emplace_back(nu + nx + 1, nx + nu);
        Llt.emplace_back(nu + nx + 1, nu);
    }

    v_AL.emplace_back(max_number_of_variables);
    v_Ggt_stripe.emplace_back(max_number_of_variables);
    v_GgLt.emplace_back(max_number_of_variables);
    v_RSQrqt_hat.emplace_back(max_number_of_variables);
    v_Llt_shift.emplace_back(max_number_of_controls);
    v_GgIt_tilde.emplace_back(info.dims.number_of_states[0]);
    v_GgLIt.emplace_back(info.dims.number_of_states[0]);
    v_HhIt.emplace_back(info.dims.number_of_states[0]);
    v_PpIt_hat.emplace_back(info.dims.number_of_states[0]);
    v_LlIt.emplace_back(info.dims.number_of_states[0]);
    v_Ggt_ineq_temp.emplace_back(max_number_of_ineq_constraints);
    v_tmp.emplace_back(max_number_of_variables);

    v_Ppt.reserve(info.dims.K);
    v_Hh.reserve(info.dims.K);
    v_RSQrqt_tilde.reserve(info.dims.K);
    v_Ggt_tilde.reserve(info.dims.K);
    v_Llt.reserve(info.dims.K);

    for (Index k = 0; k < info.dims.K; k++)
    {
        Index nu = info.dims.number_of_controls[k];
        Index nx = info.dims.number_of_states[k];
        v_Ppt.emplace_back(nx);
        v_Hh.emplace_back(nx);
        v_RSQrqt_tilde.emplace_back(nu + nx);
        v_Ggt_tilde.emplace_back(nu + nx);
        v_Llt.emplace_back(nu + nx);
    }

    PlI.emplace_back(info.dims.number_of_states[0]);
    PrI.emplace_back(info.dims.number_of_states[0]);

    Pl.reserve(info.dims.K);
    Pr.reserve(info.dims.K);

    for (Index k = 0; k < info.dims.K; k++)
    {
        Index nu = info.dims.number_of_controls[k];
        Index nx = info.dims.number_of_states[k];
        Pl.emplace_back(max_number_of_controls);
        Pr.emplace_back(max_number_of_controls);
    }

    gamma.resize(info.dims.K);
    rho.resize(info.dims.K);
};

LinsolReturnFlag OcpAugSystemSolver::solve(const ProblemInfo<OcpType> &info,
                                           Jacobian<OcpType> &jacobian, Hessian<OcpType> &hessian,
                                           const VecRealView &D_x, const VecRealView &D_s,
                                           const VecRealView &f, const VecRealView &g,
                                           VecRealView &x, VecRealView &eq_mult)
{
    MatRealView *RSQrq_hat_curr_p;
    Index rank_k;
    /////////////// recursion ///////////////
    for (Index k = info.dims.K - 1; k >= 0; --k)
    {
        const Index nu = info.dims.number_of_controls[k];
        const Index nx = info.dims.number_of_states[k];
        const Index ng = info.dims.number_of_eq_constraints[k];
        const Index ng_ineq = info.dims.number_of_ineq_constraints[k];
        const Index offset_ineq_k = info.offsets_slack[k];
        const Index offset_u = info.offsets_primal_u[k];
        const Index offset_eq_path = info.offsets_g_eq_path[k];
        const Index offset_eq_slack = info.offsets_g_eq_slack[k];
        //////// SUBSDYN
        Index gamma_k;
        if (k == info.dims.K - 1)
        {
            gamma_k = ng;
            gamma[k] = gamma_k;
            rowin(ng, 1.0, g, offset_eq_path, jacobian.Gg_eqt[k], nu + nx, 0);
            gecp(nx + nu + 1, ng, jacobian.Gg_eqt[k], 0, 0, Ggt_stripe[0], 0, 0);
            rowin(nu + nx, 1.0, f, offset_u, hessian.RSQrqt[k], nu + nx, 0);
            gecp(nx + nu + 1, nu + nx, hessian.RSQrqt[k], 0, 0, RSQrqt_tilde[k], 0, 0);
        }
        else
        {
            const Index offset_eq_dyn = info.offsets_g_eq_dyn[k];
            const Index nxp1 = info.dims.number_of_states[k + 1];
            const Index Hp1_size = gamma[k + 1] - rho[k + 1];
            if (Hp1_size + ng > nu + nx)
                return LinsolReturnFlag::NOFULL_RANK;
            gamma_k = Hp1_size + ng;
            // AL <- [BAb]^T_k P_kp1
            rowin(nxp1, 1.0, g, offset_eq_dyn, jacobian.BAbt[k], nu + nx, 0);
            gemm_nt(nu + nx + 1, nxp1, nxp1, 1.0, jacobian.BAbt[k], 0, 0, Ppt[k + 1], 0, 0, 0.0,
                    AL[0], 0, 0, AL[0], 0, 0);
            // AL[-1,:] <- AL[-1,:] + p_kp1^T
            gead(1, nxp1, 1.0, Ppt[k + 1], nxp1, 0, AL[0], nx + nu, 0);
            // RSQrqt_stripe <- AL[BA] + RSQrqt
            rowin(nu + nx, 1.0, f, offset_u, hessian.RSQrqt[k], nu + nx, 0);
            syrk_ln_mn(nu + nx + 1, nu + nx, nxp1, 1.0, AL[0], 0, 0, jacobian.BAbt[k], 0, 0, 1.0,
                       hessian.RSQrqt[k], 0, 0, RSQrqt_tilde[k], 0, 0);
            //// inequalities
            gamma[k] = gamma_k;
            // if ng[k]>0
            if (gamma_k > 0)
            {
                // if Gk nonempty
                if (ng > 0)
                {
                    // Ggt_stripe  <- Ggt_k
                    rowin(ng, 1.0, g, offset_eq_path, jacobian.Gg_eqt[k], nu + nx, 0);
                    gecp(nu + nx + 1, ng, jacobian.Gg_eqt[k], 0, 0, Ggt_stripe[0], 0, 0);
                }
                // if Hkp1 nonempty
                if (Hp1_size > 0)
                {
                    // Ggt_stripe <- [Ggt_k [BAb_k^T]H_kp1]
                    gemm_nt(nu + nx + 1, Hp1_size, nxp1, 1.0, jacobian.BAbt[k], 0, 0, Hh[k + 1], 0,
                            0, 0.0, Ggt_stripe[0], 0, ng, Ggt_stripe[0], 0, ng);
                    // Ggt_stripe[-1,ng:] <- Ggt_stripe[-1,ng:] + h_kp1^T
                    gead_transposed(1, Hp1_size, 1.0, Hh[k + 1], 0, nxp1, Ggt_stripe[0], nu + nx,
                                    ng);
                }
            }
            else
            {
                rho[k] = 0;
                rank_k = 0;
                RSQrq_hat_curr_p = &RSQrqt_tilde[k];
            }
        }
        // inequalities + inertia correction
        {
            if (ng_ineq > 0)
            {
                rowin(ng_ineq, 1.0, g, offset_eq_slack, jacobian.Gg_ineqt[k], nu + nx, 0);
                gecp(nu + nx + 1, ng_ineq, jacobian.Gg_ineqt[k], 0, 0, Ggt_ineq_temp[0], 0, 0);
                for (Index i = 0; i < ng_ineq; i++)
                {
                    Scalar scaling_factor = 1.0 / D_s(offset_ineq_k + i);
                    colsc(nu + nx + 1, scaling_factor, Ggt_ineq_temp[0], 0, i);
                }
                // add the penalty
                syrk_ln_mn(nu + nx + 1, nu + nx, ng_ineq, 1.0, Ggt_ineq_temp[0], 0, 0,
                           jacobian.Gg_ineqt[k], 0, 0, 1.0, RSQrqt_tilde[k], 0, 0, RSQrqt_tilde[k],
                           0, 0);
            }
            // inertia correction
            diaad(nu + nx, 1.0, D_x, offset_u, RSQrqt_tilde[k], 0, 0);
        }
        //////// TRANSFORM_AND_SUBSEQ
        {
            // symmetric transformation, done a little different than in paper, in order to fuse LA
            // operations LU_FACT_TRANSPOSE(Ggtstripe[:gamma_k, nu+nx+1], nu max) if(k==K-2)
            // blasfeo_print_dmat(1, gamma_k, Ggt_stripe[0], nu+nx, 0);
            lu_fact_transposed(gamma_k, nu + nx + 1, nu, rank_k, Ggt_stripe[0], Pl[k], Pr[k],
                               lu_fact_tol);

            rho[k] = rank_k;
            if (gamma_k - rank_k > 0)
            {
                // transfer eq's to next stage
                if (gamma_k - rank_k > nx)
                    return LinsolReturnFlag::NOFULL_RANK;
                getr(nx + 1, gamma_k - rank_k, Ggt_stripe[0], nu, rank_k, Hh[k], 0, 0);
            }
            if (rank_k > 0)
            {
                // Ggt_tilde_k <- Ggt_stripe[rho_k:nu+nx+1, :rho] L-T (note that this is slightly
                // different from the implementation)
                trsm_rlnn(nu - rank_k + nx + 1, rank_k, -1.0, Ggt_stripe[0], 0, 0, Ggt_stripe[0],
                          rank_k, 0, Ggt_tilde[k], 0, 0);
                // the following command copies the top block matrix (LU) to the bottom because it
                // it needed later
                gecp(rank_k, gamma_k, Ggt_stripe[0], 0, 0, Ggt_tilde[k], nu - rank_k + nx + 1, 0);
                // permutations
                trtr_l(nu + nx, RSQrqt_tilde[k], 0, 0, RSQrqt_tilde[k], 0,
                       0); // copy lower part of RSQ to upper part
                Pr[k].apply_on_rows(rank_k, &RSQrqt_tilde[k].mat()); // TODO make use of symmetry
                Pr[k].apply_on_cols(rank_k, &RSQrqt_tilde[k].mat());
                // GL <- Ggt_tilde_k @ RSQ[:rho,:nu+nx] + RSQrqt[rho:nu+nx+1, rho:] (with
                // RSQ[:rho,:nu+nx] = RSQrqt[:nu+nx,:rho]^T) GEMM_NT(nu - rank_k + nx + 1, nu + nx,
                // rank_k, 1.0, Ggt_tilde[k], 0, 0, RSQrqt_tilde[k], 0, 0, 1.0, RSQrqt_tilde_p
                // + k, rank_k, 0, GgLt[0], 0, 0); split up because valgrind was giving invalid read
                // errors when C matrix has nonzero row offset GgLt[0].print();
                gecp(nu - rank_k + nx + 1, nu + nx, RSQrqt_tilde[k], rank_k, 0, GgLt[0], 0, 0);
                gemm_nt(nu - rank_k + nx + 1, nu + nx, rank_k, 1.0, Ggt_tilde[k], 0, 0,
                        RSQrqt_tilde[k], 0, 0, 1.0, GgLt[0], 0, 0, GgLt[0], 0, 0);
                // RSQrqt_hat = GgLt[nu-rank_k + nx +1, :rank_k] * G[:rank_k, :nu+nx] +
                // GgLt[rank_k:, :]  (with G[:rank_k,:nu+nx] = Gt[:nu+nx,:rank_k]^T)
                syrk_ln_mn(nu - rank_k + nx + 1, nu + nx - rank_k, rank_k, 1.0, GgLt[0], 0, 0,
                           Ggt_tilde[k], 0, 0, 1.0, GgLt[0], 0, rank_k, RSQrqt_hat[0], 0, 0);
                // GEMM_NT(nu - rank_k + nx + 1, nu + nx - rank_k, rank_k, 1.0, GgLt[0], 0, 0,
                // Ggt_tilde[k], 0, 0, 1.0, GgLt[0], 0, rank_k, RSQrqt_hat[0], 0, 0);
                RSQrq_hat_curr_p = &RSQrqt_hat[0];
            }
            else
            {
                RSQrq_hat_curr_p = &RSQrqt_tilde[k];
            }
        }
        //////// SCHUR
        {
            if (nu - rank_k > 0)
            {
                // DLlt_k = [chol(R_hatk); Llk@chol(R_hatk)^-T]
                potrf_l_mn(nu - rank_k + nx + 1, nu - rank_k, RSQrq_hat_curr_p[0], 0, 0, Llt[k], 0,
                           0);
                if (!check_reg(nu - rank_k, &Llt[k].mat(), 0, 0))
                    return LinsolReturnFlag::INDEFINITE;
                // Pp_k = Qq_hatk - L_k^T @ Ll_k
                // SYRK_LN_MN(nx+1, nx, nu-rank_k, -1.0,Llt_p+k, nu-rank_k,0, Llt_p+k,
                // nu-rank_k,0, 1.0, RSQrq_hat_curr[0], nu-rank_k, nu-rank_k,Pp+k,0,0); // feature
                // not implmented yet
                gecp(nx + 1, nu - rank_k, Llt[k], nu - rank_k, 0, Llt_shift[0], 0,
                     0); // needless operation because feature not implemented yet
                // SYRK_LN_MN(nx + 1, nx, nu - rank_k, -1.0, Llt_shift[0], 0, 0, Llt_shift[0], 0,
                // 0, 1.0, RSQrq_hat_curr[0], nu - rank_k, nu - rank_k, Ppt[k], 0, 0);
                gecp(nx + 1, nx, RSQrq_hat_curr_p[0], nu - rank_k, nu - rank_k, Ppt[k], 0, 0);
                syrk_ln_mn(nx + 1, nx, nu - rank_k, -1.0, Llt_shift[0], 0, 0, Llt_shift[0], 0, 0,
                           1.0, Ppt[k], 0, 0, Ppt[k], 0, 0);
                // next steps are for better accuracy
                if (increased_accuracy)
                {
                    // copy eta
                    getr(nu - rank_k, gamma_k - rank_k, Ggt_stripe[0], rank_k, rank_k,
                         Ggt_stripe[0], 0, 0);
                    // blasfeo_print_dmat(gamma_k-rank_k, nu-rank_k, Ggt_stripe[0], 0,0);
                    // eta L^-T
                    trsm_rltn(gamma_k - rank_k, nu - rank_k, 1.0, Llt[k], 0, 0, Ggt_stripe[0], 0, 0,
                              Ggt_stripe[0], 0, 0);
                    // ([S^T \\ r^T] L^-T) @ (L^-1 eta^T)
                    // (eta L^-T) @ ([S^T \\ r^T] L^-T)^T
                    gemm_nt(gamma_k - rank_k, nx + 1, nu - rank_k, -1.0, Ggt_stripe[0], 0, 0,
                            Llt[k], nu - rank_k, 0, 1.0, Hh[k], 0, 0, Hh[k], 0, 0);
                    // keep (L^-1 eta^T) for forward recursion
                    getr(gamma_k - rank_k, nu - rank_k, Ggt_stripe[0], 0, 0, Ggt_tilde[k], 0,
                         rank_k);
                }
            }
            else
            {
                gecp(nx + 1, nx, RSQrq_hat_curr_p[0], 0, 0, Ppt[k], 0, 0);
            }
            trtr_l(nx, Ppt[k], 0, 0, Ppt[k], 0, 0);
        }
    }
    rankI = 0;
    //////// FIRST_STAGE
    {
        const Index nx = info.dims.number_of_states[0];
        Index gamma_I = gamma[0] - rho[0];
        if (gamma_I > nx)
        {
            return LinsolReturnFlag::NOFULL_RANK;
        }
        if (gamma_I > 0)
        {
            getr(gamma_I, nx + 1, Hh[0], 0, 0, HhIt[0], 0, 0); // transposition may be avoided
            // HhIt[0].print();
            lu_fact_transposed(gamma_I, nx + 1, nx, rankI, HhIt[0], PlI[0], PrI[0], lu_fact_tol);
            if (rankI < gamma_I)
                return LinsolReturnFlag::NOFULL_RANK;
            // PpIt_tilde <- Ggt[rankI:nx+1, :rankI] L-T (note that this is slightly different from
            // the implementation)
            trsm_rlnn(nx - rankI + 1, rankI, -1.0, HhIt[0], 0, 0, HhIt[0], rankI, 0, GgIt_tilde[0],
                      0, 0);
            // permutations
            PrI[0].apply_on_rows(rankI, &Ppt[0].mat()); // TODO make use of symmetry
            PrI[0].apply_on_cols(rankI, &Ppt[0].mat());
            // // GL <- GgIt_tilde @ Pp[:rankI,:nx] + Ppt[rankI:nx+1, rankI:] (with Pp[:rankI,:nx] =
            // Ppt[:nx,:rankI]^T) GEMM_NT(nx - rankI + 1, nx, rankI, 1.0, GgIt_tilde[0], 0, 0,
            // Ppt[0], 0, 0, 1.0, Ppt[0], rankI, 0, GgLIt[0], 0, 0); split up because valgrind was
            // giving invalid read errors when C matrix has nonzero row offset
            gecp(nx - rankI + 1, nx, Ppt[0], rankI, 0, GgLIt[0], 0, 0);
            gemm_nt(nx - rankI + 1, nx, rankI, 1.0, GgIt_tilde[0], 0, 0, Ppt[0], 0, 0, 1.0,
                    GgLIt[0], 0, 0, GgLIt[0], 0, 0);
            // // RSQrqt_hat = GgLt[nu-rank_k + nx +1, :rank_k] * G[:rank_k, :nu+nx] + GgLt[rank_k:,
            // :]  (with G[:rank_k,:nu+nx] = Gt[:nu+nx,:rank_k]^T)
            syrk_ln_mn(nx - rankI + 1, nx - rankI, rankI, 1.0, GgLIt[0], 0, 0, GgIt_tilde[0], 0, 0,
                       1.0, GgLIt[0], 0, rankI, PpIt_hat[0], 0, 0);
            // TODO skipped if nx-rankI = 0
            potrf_l_mn(nx - rankI + 1, nx - rankI, PpIt_hat[0], 0, 0, LlIt[0], 0, 0);
            if (!check_reg(nx - rankI, &LlIt[0].mat(), 0, 0))
                return LinsolReturnFlag::INDEFINITE;
        }
        else
        {
            rankI = 0;
            potrf_l_mn(nx + 1, nx, Ppt[0], 0, 0, LlIt[0], 0, 0);
            if (!check_reg(nx, &LlIt[0].mat(), 0, 0))
                return LinsolReturnFlag::INDEFINITE;
        }
    }
    ////// FORWARD_SUBSTITUTION:
    // first stage
    {
        const Index nx = info.dims.number_of_states[0];
        const Index nu = info.dims.number_of_controls[0];
        const Index offs_u = info.offsets_primal_u[0];
        const Index offs_x = info.offsets_primal_x[0];
        const Index offs_g = info.offsets_g_eq_path[0];
        // calculate xIb
        rowex(nx - rankI, -1.0, LlIt[0], nx - rankI, 0, x, offs_x + rankI);
        // assume TRSV_LTN allows aliasing, this is the case in normal BLAS
        trsv_ltn(nx - rankI, LlIt[0], 0, 0, x, offs_x + rankI, x, offs_x + rankI);
        // calculate xIa
        rowex(rankI, 1.0, GgIt_tilde[0], nx - rankI, 0, x, offs_x);
        // assume aliasing is possible for last two elements
        gemv_t(nx - rankI, rankI, 1.0, GgIt_tilde[0], 0, 0, x, offs_x + rankI, 1.0, x, offs_x, x,
               offs_x);
        //// lag
        rowex(rankI, -1.0, Ppt[0], nx, 0, eq_mult, offs_g);
        // assume aliasing is possible for last two elements
        gemv_t(nx, rankI, -1.0, Ppt[0], 0, 0, x, offs_x, 1.0, eq_mult, offs_g, eq_mult, offs_g);

        // U^-T
        trsv_lnn(rankI, HhIt[0], 0, 0, eq_mult, offs_g, eq_mult, offs_g);
        // L^-T
        trsv_unu(rankI, rankI, HhIt[0], 0, 0, eq_mult, offs_g, eq_mult, offs_g);
        PlI[0].apply_inverse(rankI, &eq_mult.vec(), offs_g);
        PrI[0].apply_inverse(rankI, &x.vec(), offs_x);
    }
    // other stages
    for (Index k = 0; k < info.dims.K; k++)
    {
        const Index nx = info.dims.number_of_states[k];
        const Index nu = info.dims.number_of_controls[k];
        const Index offs = info.offsets_primal_u[k];
        const Index offs_x = info.offsets_primal_x[k];
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
            /// calculate ukb_tilde
            // -Lkxk - lk
            rowex(numrho_k, -1.0, Llt[k], numrho_k + nx, 0, x, offs + rho_k);
            if (increased_accuracy)
            {
                gemv_n(nu - rho_k, gamma_k - rho_k, -1.0, Ggt_tilde[k], 0, rho_k, eq_mult, offs_g_k,
                       1.0, x, offs + rho_k, x, offs + rho_k);
            }
            // assume aliasing of last two eliments is allowed
            gemv_t(nx, numrho_k, -1.0, Llt[k], numrho_k, 0, x, offs_x, 1.0, x, offs + rho_k, x,
                   offs + rho_k);
            trsv_ltn(numrho_k, Llt[k], 0, 0, x, offs + rho_k, x, offs + rho_k);
        }
        /// calcualate uka_tilde
        if (rho_k > 0)
        {
            rowex(rho_k, 1.0, Ggt_tilde[k], numrho_k + nx, 0, x, offs);
            gemv_t(nx + numrho_k, rho_k, 1.0, Ggt_tilde[k], 0, 0, x, offs + rho_k, 1.0, x, offs, x,
                   offs);
            // calculate lamda_tilde_k
            // copy vk to right location
            veccp(gammamrho_k, eq_mult, offs_g_k, v_tmp[0], 0);
            veccp(gammamrho_k, v_tmp[0], 0, eq_mult, offs_g_k + rho_k);
            rowex(rho_k, -1.0, RSQrqt_tilde[k], nu + nx, 0, eq_mult, offs_g_k);
            // assume aliasing of last two eliments is allowed
            gemv_t(nu + nx, rho_k, -1.0, RSQrqt_tilde[k], 0, 0, x, offs, 1.0, eq_mult, offs_g_k,
                   eq_mult, offs_g_k);
            // nu-rank_k+nx,0
            // needless copy because feature not implemented yet in trsv_lnn
            gecp(rho_k, gamma_k, Ggt_tilde[k], nu - rho_k + nx + 1, 0, AL[0], 0, 0);
            // U^-T
            trsv_lnn(rho_k, AL[0], 0, 0, eq_mult, offs_g_k, eq_mult, offs_g_k);
            // L^-T
            trsv_unu(rho_k, gamma_k, AL[0], 0, 0, eq_mult, offs_g_k, eq_mult, offs_g_k);
            Pl[k].apply_inverse(rho_k, &eq_mult.vec(), offs_g_k);
            Pr[k].apply_inverse(rho_k, &x.vec(), offs);
        }
        if (ng_ineq > 0)
        {
            gemv_t(nu + nx, ng_ineq, 1.0, jacobian.Gg_ineqt[k], 0, 0, x, offs, 1.0, g, offs_eq_ineq,
                   eq_mult, offs_eq_ineq);
            eq_mult.block(ng_ineq, offs_eq_ineq) =
                eq_mult.block(ng_ineq, offs_eq_ineq) / D_s.block(ng_ineq, offs_slack);
        }
        if (k != info.dims.K - 1)
        {
            const Index offs_dyn_eq_k = info.offsets_g_eq_dyn[k];
            const Index nxp1 = info.dims.number_of_states[k + 1];
            const Index nup1 = info.dims.number_of_controls[k + 1];
            const Index offsp1 = info.offsets_primal_u[k + 1];
            const Index offsxp1 = info.offsets_primal_x[k + 1];
            const Index offs_g_kp1 = info.offsets_g_eq_path[k + 1];
            const Index gammamrho_kp1 = gamma[k + 1] - rho[k + 1];
            // calculate xkp1
            rowex(nxp1, 1.0, jacobian.BAbt[k], nu + nx, 0, x, offsxp1);
            gemv_t(nu + nx, nxp1, 1.0, jacobian.BAbt[k], 0, 0, x, offs, 1.0, x, offsxp1, x,
                   offsxp1);
            // calculate lam_dyn xp1
            rowex(nxp1, 1.0, Ppt[k + 1], nxp1, 0, eq_mult, offs_dyn_eq_k);
            gemv_t(nxp1, nxp1, 1.0, Ppt[k + 1], 0, 0, x, offsxp1, 1.0, eq_mult, offs_dyn_eq_k,
                   eq_mult, offs_dyn_eq_k);
            gemv_t(gammamrho_kp1, nxp1, 1.0, Hh[k + 1], 0, 0, eq_mult, offs_g_kp1, 1.0, eq_mult,
                   offs_dyn_eq_k, eq_mult, offs_dyn_eq_k);
        }
    }
    return LinsolReturnFlag::SUCCESS;
}
LinsolReturnFlag OcpAugSystemSolver::solve(const ProblemInfo<OcpType> &info,
                                           Jacobian<OcpType> &jacobian, Hessian<OcpType> &hessian,
                                           const VecRealView &D_x, const VecRealView &D_eq,
                                           const VecRealView &D_s, const VecRealView &f,
                                           const VecRealView &g, VecRealView &x,
                                           VecRealView &eq_mult)
{
    MatRealView *RSQrq_hat_curr_p;
    for (Index k = info.dims.K - 1; k >= 0; --k)
    {
        const Index nu = info.dims.number_of_controls[k];
        const Index nx = info.dims.number_of_states[k];
        const Index ng = info.dims.number_of_eq_constraints[k];
        const Index ng_ineq = info.dims.number_of_ineq_constraints[k];
        const Index offs_ineq_k = info.offsets_slack[k];
        const Index offset_u = info.offsets_primal_u[k];
        const Index offset_eq_k = info.offsets_eq[k];
        const Index offset_g_eq_k = info.offsets_g_eq_path[k];
        const Index offset_g_ineq_k = info.offsets_g_eq_slack[k];
        // const fatrop_int offs_g_ineq_k = offs_g_ineq_p[k];
        //////// SUBSDYN
        if (k == info.dims.K - 1)
        {
            rowin(nu + nx, 1.0, f, offset_u, hessian.RSQrqt[k], nu + nx, 0);
            gecp(nx + nu + 1, nu + nx, hessian.RSQrqt[k], 0, 0, RSQrqt_tilde[k], 0, 0);
        }
        else
        {
            const Index offset_eq_dyn = info.offsets_g_eq_dyn[k];
            const Index nxp1 = info.dims.number_of_states[k + 1];
            // AL <- [BAb]^T_k P_kp1
            rowin(nxp1, 1.0, g, offset_eq_dyn, jacobian.BAbt[k], nu + nx, 0);
            gemm_nt(nu + nx + 1, nxp1, nxp1, 1.0, jacobian.BAbt[k], 0, 0, Ppt[k + 1], 0, 0, 0.0,
                    AL[0], 0, 0, AL[0], 0, 0);
            // AL[-1,:] <- AL[-1,:] + p_kp1^T
            gead(1, nxp1, 1.0, Ppt[k + 1], nxp1, 0, AL[0], nx + nu, 0);
            // RSQrqt_stripe <- AL[BA] + RSQrqt
            rowin(nu + nx, 1.0, f, offset_u, hessian.RSQrqt[k], nu + nx, 0);
            syrk_ln_mn(nu + nx + 1, nu + nx, nxp1, 1.0, AL[0], 0, 0, jacobian.BAbt[k], 0, 0, 1.0,
                       hessian.RSQrqt[k], 0, 0, RSQrqt_tilde[k], 0, 0);
        }
        // equality penalty
        {
            rowin(ng, 1.0, g, offset_g_eq_k, jacobian.Gg_eqt[k], nu + nx, 0);
            gecp(nu + nx + 1, ng, jacobian.Gg_eqt[k], 0, 0, Ggt_stripe[0], 0, 0);
            for (Index i = 0; i < ng; i++)
            {
                Scalar scaling_factor = 1.0 / D_eq(offset_eq_k + i);
                colsc(nu + nx + 1, scaling_factor, Ggt_stripe[0], 0, i);
            }
            // add the penalty
            syrk_ln_mn(nu + nx + 1, nu + nx, ng, 1.0, Ggt_stripe[0], 0, 0, jacobian.Gg_eqt[k], 0, 0,
                       1.0, RSQrqt_tilde[k], 0, 0, RSQrqt_tilde[k], 0, 0);
        }
        // inequalities + inertia correction
        {
            if (ng_ineq > 0)
            {
                rowin(ng_ineq, 1.0, g, offset_g_ineq_k, jacobian.Gg_ineqt[k], nu + nx, 0);
                gecp(nu + nx + 1, ng_ineq, jacobian.Gg_ineqt[k], 0, 0, Ggt_ineq_temp[0], 0, 0);
                for (Index i = 0; i < ng_ineq; i++)
                {
                    Scalar scaling_factor = 1.0 / D_s(offs_ineq_k + i);
                    colsc(nu + nx + 1, scaling_factor, Ggt_ineq_temp[0], 0, i);
                }
                // add the penalty
                syrk_ln_mn(nu + nx + 1, nu + nx, ng_ineq, 1.0, Ggt_ineq_temp[0], 0, 0,
                           jacobian.Gg_ineqt[k], 0, 0, 1.0, RSQrqt_tilde[k], 0, 0, RSQrqt_tilde[k],
                           0, 0);
            }
            // inertia correction
            diaad(nu + nx, 1.0, D_x, offset_u, RSQrqt_tilde[k], 0, 0);
        }

        //////// TRANSFORM_AND_SUBSEQ
        {
            RSQrq_hat_curr_p = &RSQrqt_tilde[k];
        }
        //////// SCHUR
        {
            // DLlt_k = [chol(R_hatk); Llk@chol(R_hatk)^-T]
            potrf_l_mn(nu + nx + 1, nu, *RSQrq_hat_curr_p, 0, 0, Llt[k], 0, 0);
            if (!check_reg(nu, &Llt[k].mat(), 0, 0))
                return LinsolReturnFlag::INDEFINITE;
            // Pp_k = Qq_hatk - L_k^T @ Ll_k
            // SYRK_LN_MN(nx+1, nx, nu-rank_k, -1.0,Llt_p+k, nu-rank_k,0, Llt_p+k, nu-rank_k,0, 1.0,
            // RSQrq_hat_curr_p, nu-rank_k, nu-rank_k,Pp+k,0,0); // feature not implmented yet
            gecp(nx + 1, nu, Llt[k], nu, 0, Llt_shift[0], 0,
                 0); // needless operation because feature not implemented yet
            syrk_ln_mn(nx + 1, nx, nu, -1.0, Llt_shift[0], 0, 0, Llt_shift[0], 0, 0, 1.0,
                       *RSQrq_hat_curr_p, nu, nu, Ppt[k], 0, 0);
        }
        trtr_l(nx, Ppt[k], 0, 0, Ppt[k], 0, 0);
    }
    //////// FIRST_STAGE
    {
        const Index nx = info.dims.number_of_states[0];
        {
            potrf_l_mn(nx + 1, nx, Ppt[0], 0, 0, LlIt[0], 0, 0);
            if (!check_reg(nx, &LlIt[0].mat(), 0, 0))
                return LinsolReturnFlag::INDEFINITE;
        }
    }
    ////// FORWARD_SUBSTITUTION:
    // first stage
    {
        const Index nx = info.dims.number_of_states[0];
        const Index nu = info.dims.number_of_controls[0];
        const Index offs_x = info.offsets_primal_x[0];
        // calculate xIb
        rowex(nx, -1.0, LlIt[0], nx, 0, x, offs_x);
        // assume TRSV_LTN allows aliasing, this is the case in normal BLAS
        trsv_ltn(nx, LlIt[0], 0, 0, x, offs_x, x, offs_x);
    }
    for (Index k = 0; k < info.dims.K; k++)
    {
        const Index nx = info.dims.number_of_states[k];
        const Index nu = info.dims.number_of_controls[k];
        const Index offs = info.offsets_primal_u[k];
        const Index offs_x = info.offsets_primal_x[k];
        rowex(nu, -1.0, Llt[k], nu + nx, 0, x, offs);
        gemv_t(nx, nu, -1.0, Llt[k], nu, 0, x, offs_x, 1.0, x, offs, x, offs);
        trsv_ltn(nu, Llt[k], 0, 0, x, offs, x, offs);
        if (k != info.dims.K - 1)
        {
            const Index nxp1 = info.dims.number_of_states[k + 1];
            const Index nup1 = info.dims.number_of_controls[k + 1];
            const Index offs_x_p1 = info.offsets_primal_x[k + 1];
            const Index offs_dyn_eq_k = info.offsets_g_eq_dyn[k];
            // calculate xkp1
            rowex(nxp1, 1.0, jacobian.BAbt[k], nu + nx, 0, x, offs_x_p1);
            gemv_t(nu + nx, nxp1, 1.0, jacobian.BAbt[k], 0, 0, x, offs, 1.0, x, offs_x_p1, x,
                   offs_x_p1);
            // calculate lam_dyn xp1
            rowex(nxp1, 1.0, Ppt[k + 1], nxp1, 0, eq_mult, offs_dyn_eq_k);
            gemv_t(nxp1, nxp1, 1.0, Ppt[k + 1], 0, 0, x, offs_x_p1, 1.0, eq_mult, offs_dyn_eq_k,
                   eq_mult, offs_dyn_eq_k);
        }
        const Index ng = info.dims.number_of_eq_constraints[k];
        const Index offs_g_eq_k = info.offsets_g_eq_path[k];
        const Index offs_eq_k = info.offsets_eq[k];
        if (ng > 0)
        {
            gemv_t(nu + nx, ng, 1.0, jacobian.Gg_eqt[k], 0, 0, x, offs, 1.0, g, offs_g_eq_k,
                   eq_mult, offs_g_eq_k);
            eq_mult.block(ng, offs_g_eq_k) =
                eq_mult.block(ng, offs_g_eq_k) / D_eq.block(ng, offs_eq_k);
        }
        const Index ng_ineq = info.dims.number_of_ineq_constraints[k];
        const Index offs_slack = info.offsets_slack[k];
        const Index offs_eq_ineq = info.offsets_g_eq_slack[k];
        if (ng_ineq > 0)
        {
            gemv_t(nu + nx, ng_ineq, 1.0, jacobian.Gg_ineqt[k], 0, 0, x, offs, 1.0, g, offs_eq_ineq,
                   eq_mult, offs_eq_ineq);
            eq_mult.block(ng_ineq, offs_eq_ineq) =
                eq_mult.block(ng_ineq, offs_eq_ineq) / D_s.block(ng_ineq, offs_slack);
        }
    }
    return LinsolReturnFlag::SUCCESS;
}

LinsolReturnFlag OcpAugSystemSolver::solve_rhs(const ProblemInfo<OcpType> &info,
                                               const Jacobian<OcpType> &jacobian,
                                               const Hessian<OcpType> &hessian,
                                               const VecRealView &D_s, const VecRealView &f,
                                               const VecRealView &g, VecRealView &x,
                                               VecRealView &eq_mult)
{
    VecRealView *v_RSQrq_hat_curr_p;
    Index rank_k;
    /////////////// recursion ///////////////

    for (Index k = info.dims.K - 1; k >= 0; --k)
    {
        const Index nu = info.dims.number_of_controls[k];
        const Index nx = info.dims.number_of_states[k];
        const Index ng = info.dims.number_of_eq_constraints[k];
        const Index ng_ineq = info.dims.number_of_ineq_constraints[k];
        const Index offset_ineq_k = info.offsets_slack[k];
        const Index offs_g_ineq_k = info.offsets_g_eq_slack[k];
        const Index offs_g_k = info.offsets_g_eq_path[k];
        const Index offs = info.offsets_primal_u[k];
        //         //////// SUBSDYN
        Index gamma_k;
        if (k == info.dims.K - 1)
        {
            gamma_k = ng;
            gamma[k] = gamma_k;
            veccp(ng, g, offs_g_k, v_Ggt_stripe[0], 0);
            veccp(nu + nx, f, offs, v_RSQrqt_tilde[k], 0);
        }
        else
        {
            const Index offs_dyn_k = info.offsets_g_eq_dyn[k];
            const Index nxp1 = info.dims.number_of_states[k + 1];
            const Index Hp1_size = gamma[k + 1] - rho[k + 1];
            gamma_k = Hp1_size + ng;
            gemv_n(nxp1, nxp1, 1.0, Ppt[k + 1], 0, 0, g, offs_dyn_k, 0.0, v_AL[0], 0, v_AL[0], 0);
            axpy(nxp1, 1.0, v_Ppt[k + 1], 0, v_AL[0], 0, v_AL[0], 0);
            gemv_n(nu + nx, nxp1, 1.0, jacobian.BAbt[k], 0, 0, v_AL[0], 0, 1.0, f, offs,
                   v_RSQrqt_tilde[k], 0);
            if (gamma_k > 0)
            {
                if (ng > 0)
                {
                    veccp(ng, g, offs_g_k, v_Ggt_stripe[0], 0);
                }
                if (Hp1_size > 0)
                {
                    gemv_n(Hp1_size, nxp1, 1.0, Hh[k + 1], 0, 0, g, offs_dyn_k, 0.0,
                           v_Ggt_stripe[0], ng, v_Ggt_stripe[0], ng);
                    axpy(Hp1_size, 1.0, v_Hh[k + 1], 0, v_Ggt_stripe[0], ng, v_Ggt_stripe[0], ng);
                }
            }
            else
            {
                rank_k = 0;
                v_RSQrq_hat_curr_p = &v_RSQrqt_tilde[k];
            }
        }
        if (ng_ineq > 0)
        {
            for (Index i = 0; i < ng_ineq; i++)
            {
                Scalar scaling_factor = D_s(offset_ineq_k + i);
                Scalar grad_barrier = g(offs_g_ineq_k + i);
                v_Ggt_ineq_temp[0](i) = grad_barrier / scaling_factor;
            }
            gemv_n(nu + nx, ng_ineq, 1.0, jacobian.Gg_ineqt[k], 0, 0, v_Ggt_ineq_temp[0], 0, 1.0,
                   v_RSQrqt_tilde[k], 0, v_RSQrqt_tilde[k], 0);
        }
        {
            rank_k = rho[k];
            gecp(rank_k, gamma_k, Ggt_tilde[k], nu - rank_k + nx + 1, 0, Ggt_stripe[0], 0, 0);
            Pl[k].apply(rank_k, &v_Ggt_stripe[0].vec(), 0);
            trsv_utu(rank_k, Ggt_stripe[0], 0, 0, v_Ggt_stripe[0], 0, v_Ggt_stripe[0], 0);
            gemv_t(rank_k, gamma_k - rank_k, -1.0, Ggt_stripe[0], 0, rank_k, v_Ggt_stripe[0], 0,
                   1.0, v_Ggt_stripe[0], rank_k, v_Ggt_stripe[0], rank_k);

            if (gamma_k - rank_k > 0)
            {
                veccp(gamma_k - rank_k, v_Ggt_stripe[0], rank_k, v_Hh[k], 0);
            }
            if (rank_k > 0)
            {
                veccpsc(rank_k, -1.0, v_Ggt_stripe[0], 0, v_Ggt_tilde[k], 0);
                trsv_ltn(rank_k, Ggt_stripe[0], 0, 0, v_Ggt_tilde[k], 0, v_Ggt_tilde[k], 0);
                Pr[k].apply(rank_k, &v_RSQrqt_tilde[k].vec(), 0);
                veccp(nu + nx, v_RSQrqt_tilde[k], 0, v_GgLt[0], 0);
                gemv_n(nu + nx, rank_k, 1.0, RSQrqt_tilde[k], 0, 0, v_Ggt_tilde[k], 0, 1.0,
                       v_GgLt[0], 0, v_GgLt[0], 0);
                gemv_n(nu + nx - rank_k, rank_k, 1.0, Ggt_tilde[k], 0, 0, v_GgLt[0], 0, 1.0,
                       v_GgLt[0], rank_k, v_RSQrqt_hat[0], 0);
                v_RSQrq_hat_curr_p = &v_RSQrqt_hat[0];
            }
            else
            {
                v_RSQrq_hat_curr_p = &v_RSQrqt_tilde[k];
            }
        }
        //         //////// SCHUR
        {
            if (nu - rank_k > 0)
            {
                trsv_lnn(nu - rank_k, Llt[k], 0, 0, *v_RSQrq_hat_curr_p, 0, v_Llt[k], 0);
                gecp(nx + 1, nu - rank_k, Llt[k], nu - rank_k, 0, Llt_shift[0], 0, 0);
                veccp(nu - rank_k, v_Llt[k], 0, v_Llt_shift[0], 0);
                veccp(nx, *v_RSQrq_hat_curr_p, nu - rank_k, v_Ppt[k], 0);
                gemv_n(nx, nu - rank_k, -1.0, Llt_shift[0], 0, 0, v_Llt_shift[0], 0, 1.0, v_Ppt[k],
                       0, v_Ppt[k], 0);
                if (increased_accuracy)
                {
                    gemv_t(nu - rank_k, gamma_k - rank_k, -1.0, Ggt_tilde[k], 0, rank_k, v_Llt[k],
                           0, 1.0, v_Hh[k], 0, v_Hh[k], 0);
                }
            }
            else
            {
                veccp(nx, *v_RSQrq_hat_curr_p, 0, v_Ppt[k], 0);
            }
        }
    }
    {
        const Index nx = info.dims.number_of_states[0];
        Index gamma_I = gamma[0] - rho[0];
        if (gamma_I > 0)
        {
            veccp(gamma_I, v_Hh[0], 0, v_HhIt[0], 0);
            PlI[0].apply(rankI, &v_HhIt[0].vec(), 0);
            trsv_utu(rankI, HhIt[0], 0, 0, v_HhIt[0], 0, v_HhIt[0], 0);
            gemv_t(rankI, gamma_I - rankI, -1.0, HhIt[0], 0, rankI, v_HhIt[0], 0, 1.0, v_HhIt[0],
                   rankI, v_HhIt[0], rankI);
            veccpsc(rankI, -1.0, v_HhIt[0], 0, v_GgIt_tilde[0], 0);
            trsv_ltn(rankI, HhIt[0], 0, 0, v_GgIt_tilde[0], 0, v_GgIt_tilde[0], 0);
            PrI[0].apply(rankI, &v_Ppt[0].vec(), 0);
            veccp(nx, v_Ppt[0], 0, v_GgLIt[0], 0);
            gemv_n(nx, rankI, 1.0, Ppt[0], 0, 0, v_GgIt_tilde[0], 0, 1.0, v_GgLIt[0], 0, v_GgLIt[0],
                   0);
            gemv_n(nx - rankI, rankI, 1.0, GgIt_tilde[0], 0, 0, v_GgLIt[0], 0, 1.0, v_GgLIt[0],
                   rankI, v_PpIt_hat[0], 0);
            trsv_lnn(nx - rankI, LlIt[0], 0, 0, v_PpIt_hat[0], 0, v_LlIt[0], 0);
        }
        else
        {
            trsv_lnn(nx, LlIt[0], 0, 0, v_Ppt[0], 0, v_LlIt[0], 0);
        }
    }
    {
        const Index nx = info.dims.number_of_states[0];
        const Index nu = info.dims.number_of_controls[0];
        const Index offs_u = info.offsets_primal_u[0];
        const Index offs_x = info.offsets_primal_x[0];
        const Index offs_g = info.offsets_g_eq_path[0];
        veccpsc(nx - rankI, -1.0, v_LlIt[0], 0, x, offs_x + rankI);
        trsv_ltn(nx - rankI, LlIt[0], 0, 0, x, offs_x + rankI, x, offs_x + rankI);
        veccp(rankI, v_GgIt_tilde[0], 0, x, offs_x);
        gemv_t(nx - rankI, rankI, 1.0, GgIt_tilde[0], 0, 0, x, offs_x + rankI, 1.0, x, offs_x, x,
               offs_x);
        veccpsc(rankI, -1.0, v_Ppt[0], 0, eq_mult, offs_g);
        gemv_t(nx, rankI, -1.0, Ppt[0], 0, 0, x, nu, 1.0, eq_mult, offs_g, eq_mult, offs_g);
        trsv_lnn(rankI, HhIt[0], 0, 0, eq_mult, offs_g, eq_mult, offs_g);
        trsv_unu(rankI, rankI, HhIt[0], 0, 0, eq_mult, offs_g, eq_mult, offs_g);
        PlI[0].apply_inverse(rankI, &eq_mult.vec(), offs_g);
        PrI[0].apply_inverse(rankI, &x.vec(), offs_x);
    }
    for (Index k = 0; k < info.dims.K; k++)
    {

        const Index nx = info.dims.number_of_states[k];
        const Index nu = info.dims.number_of_controls[k];
        const Index offs = info.offsets_primal_u[k];
        const Index offs_x = info.offsets_primal_x[k];
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
            veccpsc(numrho_k, -1.0, v_Llt[k], 0, x, offs + rho_k);
            if (increased_accuracy)
            {
                gemv_n(nu - rho_k, gamma_k - rho_k, -1.0, Ggt_tilde[k], 0, rho_k, eq_mult, offs_g_k,
                       1.0, x, offs + rho_k, x, offs + rho_k);
            }
            gemv_t(nx, numrho_k, -1.0, Llt[k], numrho_k, 0, x, offs_x, 1.0, x, offs + rho_k, x,
                   offs + rho_k);
            trsv_ltn(numrho_k, Llt[k], 0, 0, x, offs + rho_k, x, offs + rho_k);
        }
        //         /// calcualate uka_tilde
        if (rho_k > 0)
        {
            // ROWEX(rho_k, 1.0, Ggt_tilde[k], numrho_k + nx, 0, ux[0], offs);
            veccp(rho_k, v_Ggt_tilde[k], 0, x, offs);
            gemv_t(nx + numrho_k, rho_k, 1.0, Ggt_tilde[k], 0, 0, x, offs + rho_k, 1.0, x, offs, x,
                   offs);
            veccp(gammamrho_k, eq_mult, offs_g_k, v_tmp[0], 0);
            veccp(gammamrho_k, v_tmp[0], 0, eq_mult, offs_g_k + rho_k);
            veccpsc(rho_k, -1.0, v_RSQrqt_tilde[k], 0, eq_mult, offs_g_k);
            gemv_t(nu + nx, rho_k, -1.0, RSQrqt_tilde[k], 0, 0, x, offs, 1.0, eq_mult, offs_g_k,
                   eq_mult, offs_g_k);
            gecp(rho_k, gamma_k, Ggt_tilde[k], nu - rho_k + nx + 1, 0, AL[0], 0, 0);
            trsv_lnn(rho_k, AL[0], 0, 0, eq_mult, offs_g_k, eq_mult, offs_g_k);
            trsv_unu(rho_k, gamma_k, AL[0], 0, 0, eq_mult, offs_g_k, eq_mult, offs_g_k);
            Pl[k].apply_inverse(rho_k, &eq_mult.vec(), offs_g_k);
            Pr[k].apply_inverse(rho_k, &x.vec(), offs);
        }
        if (ng_ineq > 0)
        {
            gemv_t(nu + nx, ng_ineq, 1.0, jacobian.Gg_ineqt[k], 0, 0, x, offs, 1.0, g, offs_eq_ineq,
                   eq_mult, offs_eq_ineq);
            eq_mult.block(ng_ineq, offs_eq_ineq) =
                eq_mult.block(ng_ineq, offs_eq_ineq) / D_s.block(ng_ineq, offs_slack);
        }
        if (k != info.dims.K - 1)
        {
            const Index nxp1 = info.dims.number_of_states[k + 1];
            const Index nup1 = info.dims.number_of_controls[k + 1];
            const Index offsp1 = info.offsets_primal_u[k + 1];
            const Index offsxp1 = info.offsets_primal_x[k + 1];
            const Index offs_g_kp1 = info.offsets_g_eq_path[k + 1];
            const Index offs_dyn_k = info.offsets_g_eq_dyn[k];
            const Index gammamrho_kp1 = gamma[k + 1] - rho[k + 1];
            veccp(nxp1, g, offs_dyn_k, x, offsxp1);
            gemv_t(nu + nx, nxp1, 1.0, jacobian.BAbt[k], 0, 0, x, offs, 1.0, x, offsxp1, x,
                   offsxp1);
            veccp(nxp1, v_Ppt[k + 1], 0, eq_mult, offs_dyn_k);
            gemv_t(nxp1, nxp1, 1.0, Ppt[k + 1], 0, 0, x, offsxp1, 1.0, eq_mult, offs_dyn_k, eq_mult,
                   offs_dyn_k);
            gemv_t(gammamrho_kp1, nxp1, 1.0, Hh[k + 1], 0, 0, eq_mult, offs_g_kp1, 1.0, eq_mult,
                   offs_dyn_k, eq_mult, offs_dyn_k);
        }
    }
    return LinsolReturnFlag::SUCCESS;
}
LinsolReturnFlag OcpAugSystemSolver::solve_rhs(const ProblemInfo<OcpType> &info,
                                               const Jacobian<OcpType> &jacobian,
                                               const Hessian<OcpType> &hessian,
                                               const VecRealView &D_eq, const VecRealView &D_s,
                                               const VecRealView &f, const VecRealView &g,
                                               VecRealView &x, VecRealView &eq_mult)
{
    VecRealView *v_RSQrq_hat_curr_p;
    for (Index k = info.dims.K - 1; k >= 0; --k)
    {
        const Index offs_ux_k = info.offsets_primal_u[k];
        const Index nu = info.dims.number_of_controls[k];
        const Index nx = info.dims.number_of_states[k];
        const Index ng = info.dims.number_of_eq_constraints[k];
        const Index ng_ineq = info.dims.number_of_ineq_constraints[k];
        const Index offs_g_dyn = info.offsets_g_eq_dyn[k];
        const Index offs_g_eq = info.offsets_g_eq_path[k];
        const Index offs_ge_eq_ineq = info.offsets_g_eq_slack[k];
        //     //////// SUBSDYN
        if (k == info.dims.K - 1)
        {
            veccp(nu + nx, f, offs_ux_k, v_RSQrqt_tilde[k], 0);
        }
        else
        {
            const Index nxp1 = info.dims.number_of_states[k + 1];
            gemv_n(nxp1, nxp1, 1.0, Ppt[k + 1], 0, 0, g, offs_g_dyn, 0.0, v_AL[0], 0, v_AL[0], 0);
            axpy(nxp1, 1.0, v_Ppt[k + 1], 0, v_AL[0], 0, v_AL[0], 0);
            gemv_n(nu + nx, nxp1, 1.0, jacobian.BAbt[k], 0, 0, v_AL[0], 0, 1.0, f, offs_ux_k,
                   v_RSQrqt_tilde[k], 0);
        }
        if (ng > 0)
        {
            const Index offs_eq_k = info.offsets_eq[k];
            for (Index i = 0; i < ng; i++)
            {
                Scalar scaling_factor = D_eq(offs_eq_k + i);
                v_Ggt_stripe[0](i) = g(offs_g_eq + i) / scaling_factor;
            }
            gemv_n(nu + nx, ng, 1.0, jacobian.Gg_eqt[k], 0, 0, v_Ggt_stripe[0], 0, 1.0,
                   v_RSQrqt_tilde[k], 0, v_RSQrqt_tilde[k], 0);
        }
        if (ng_ineq > 0)
        {
            const Index offs_ineq_k = info.offsets_slack[k];
            for (Index i = 0; i < ng_ineq; i++)
            {
                Scalar scaling_factor = D_s(offs_ineq_k + i);
                v_Ggt_ineq_temp[0](i) = g(offs_ge_eq_ineq + i) / scaling_factor;
            }
            gemv_n(nu + nx, ng_ineq, 1.0, jacobian.Gg_ineqt[k], 0, 0, v_Ggt_ineq_temp[0], 0, 1.0,
                   v_RSQrqt_tilde[k], 0, v_RSQrqt_tilde[k], 0);
        }
        {
            v_RSQrq_hat_curr_p = &v_RSQrqt_tilde[k];
        }
        {
            trsv_lnn(nu, Llt[k], 0, 0, *v_RSQrq_hat_curr_p, 0, v_Llt[k], 0);
            veccp(nu, v_Llt[k], 0, v_Llt_shift[0], 0);
            gemv_n(nx, nu, -1.0, Llt[k], nu, 0, v_Llt_shift[0], 0, 1.0, v_RSQrqt_tilde[k], nu,
                   v_Ppt[k], 0);
        }
    }
    {
        const Index nx = info.dims.number_of_states[0];
        {
            trsv_lnn(nx, LlIt[0], 0, 0, v_Ppt[0], 0, v_LlIt[0], 0);
        }
    }
    {
        const Index nx = info.dims.number_of_states[0];
        const Index nu = info.dims.number_of_controls[0];
        const Index offs_x = info.offsets_primal_x[0];
        veccpsc(nx, -1.0, v_LlIt[0], 0, x, offs_x);
        trsv_ltn(nx, LlIt[0], 0, 0, x, offs_x, x, offs_x);
    }
    for (Index k = 0; k < info.dims.K; k++)
    {
        const Index nx = info.dims.number_of_states[k];
        const Index nu = info.dims.number_of_controls[k];
        const Index offs = info.offsets_primal_u[k];
        const Index offs_x = info.offsets_primal_x[k];
        const Index offs_dyn_eq_k = info.offsets_g_eq_dyn[k];
        veccpsc(nu, -1.0, v_Llt[k], 0, x, offs);
        gemv_t(nx, nu, -1.0, Llt[k], nu, 0, x, offs_x, 1.0, x, offs, x, offs);
        trsv_ltn(nu, Llt[k], 0, 0, x, offs, x, offs);
        if (k != info.dims.K - 1)
        {
            const Index nxp1 = info.dims.number_of_states[k + 1];
            const Index offsp1 = info.offsets_primal_u[k + 1];
            const Index offs_x_p1 = info.offsets_primal_x[k + 1];
            veccp(nxp1, g, offs_dyn_eq_k, x, offs_x_p1);
            gemv_t(nu + nx, nxp1, 1.0, jacobian.BAbt[k], 0, 0, x, offs, 1.0, x, offs_x_p1, x,
                   offs_x_p1);
            veccp(nxp1, v_Ppt[k + 1], 0, eq_mult, offs_dyn_eq_k);
            gemv_t(nxp1, nxp1, 1.0, Ppt[k + 1], 0, 0, x, offs_x_p1, 1.0, eq_mult,
                   offs_dyn_eq_k, eq_mult, offs_dyn_eq_k);
        }
    }
    // // calculate lam_eq xk
    for (Index k = 0; k < info.dims.K; k++)
    {
        const Index nx = info.dims.number_of_states[k];
        const Index nu = info.dims.number_of_controls[k];
        const Index ng = info.dims.number_of_eq_constraints[k];
        const Index offs = info.offsets_primal_u[k];
        const Index offs_g_k = info.offsets_g_eq_path[k];
        const Index offs_eq = info.offsets_eq[k];
        if (ng > 0)
        {
            gemv_t(nu + nx, ng, 1.0, jacobian.Gg_eqt[k], 0, 0, x, offs, 1.0, g, offs_g_k,
                   eq_mult, offs_g_k);
            eq_mult.block(ng, offs_g_k) =
                eq_mult.block(ng, offs_g_k) / D_eq.block(ng, offs_eq);
        }
    }

    for (Index k = 0; k < info.dims.K; k++)
    {
        const Index nx = info.dims.number_of_states[k];
        const Index nu = info.dims.number_of_controls[k];
        const Index ng_ineq = info.dims.number_of_ineq_constraints[k];
        const Index offs = info.offsets_primal_u[k];
        const Index offs_gineq_k = info.offsets_g_eq_slack[k];
        const Index offs_slack = info.offsets_slack[k];
        if (ng_ineq > 0)
        {
            gemv_t(nu + nx, ng_ineq, 1.0, jacobian.Gg_ineqt[k], 0, 0, x, offs, 1.0, g, offs_gineq_k,
                   eq_mult, offs_gineq_k);
            eq_mult.block(ng_ineq, offs_gineq_k) =
                eq_mult.block(ng_ineq, offs_gineq_k) / D_s.block(ng_ineq, offs_slack);
        }
    }
    return LinsolReturnFlag::SUCCESS;
}