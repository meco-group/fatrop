#include "LineSearchDDP.hpp"
using namespace fatrop;
int LineSearchDDP::TryStep(double alpha_primal, double alpha_dual) const
{
    int K = OCP_->K;
    // make variables local for efficiency
    MAT *BAbt_p = (MAT *)OCP_->BAbt;
    MAT *Ggt_ineq_p = (MAT *)OCP_->Ggt_ineq;
    MAT *Ppt_p = (MAT *) ocplsriccati_->Ppt;
    MAT * Hh_p = (MAT *) ocplsriccati_->Hh;
    MAT * AL_p = (MAT *) ocplsriccati_->AL;
    MAT * RSQrqt_tilde_p = (MAT *) ocplsriccati_->RSQrqt_tilde;
    MAT * Ggt_tilde_p = (MAT *) ocplsriccati_->Ggt_tilde;
    PMAT * Pl_p = (PMAT *) ocplsriccati_->Pl;
    PMAT * Pr_p = (PMAT *) ocplsriccati_->Pr;
    MAT * Llt_p = (MAT *) ocplsriccati_->Llt;
    MAT * GgIt_tilde_p = (MAT *) ocplsriccati_->GgIt_tilde;
    MAT * HhIt_p = (MAT *) ocplsriccati_->HhIt;
    MAT * LlIt_p = (MAT *) ocplsriccati_->LlIt;
    PMAT * PlI_p = (PMAT *) ocplsriccati_->PlI;
    PMAT * PrI_p = (PMAT *) ocplsriccati_->PrI;
    VEC * delta_ux_p = (VEC *) fatropdata_->delta_x;
    VEC * next_ux_p = (VEC *) fatropdata_->x_next;
    VEC * curr_ux_p = (VEC *) fatropdata_->x_curr;
    // VEC * lam_p = (VEC *) fatropdata_->lam_calc;
    VEC * lam_curr_p = (VEC *) fatropdata_->lam_curr;
    VEC * s_p = (VEC *) fatropdata_->s_curr;
    VEC * zL_p = (VEC *) fatropdata_->zL_curr;
    VEC * zU_p = (VEC *) fatropdata_->zU_curr;
    // VEC * delta_zL_p = (VEC *) fatropdata_->delta_zL;
    // VEC * delta_zU_p = (VEC *) fatropdata_->delta_zU;
    VEC * lower_p = (VEC *) fatropdata_->s_lower;
    VEC * upper_p = (VEC *) fatropdata_->s_upper;
    VEC * delta_s_p = (VEC *) fatropdata_->delta_s;
    int * nu_p = (int *) OCP_->nu;
    int * nx_p = (int *) OCP_->nx;
    int * ng_ineq_p = (int *) OCP_->ng_ineq;
    int * gamma_p = (int *) ocplsriccati_->gamma;
    int * rho_p = (int *) ocplsriccati_->rho;
    int *offs_ineq_p = (int *)OCP_->aux.ineq_offs.data();
    int *offs_g_ineq_p = (int *)OCP_->aux.g_ineq_offs.data();
    int rankI = ocplsriccati_->lastused_.rankI;
    double inertia_correction = ocplsriccati_->lastused_.inertia_correction;
    double kappa_d = ocplsriccati_->lastused_.kappa_d;
    double mu = ocplsriccati_->lastused_.mu;
    // perform a forward DDP pass
    ////// FORWARD_SUBSTITUTION:
    // first stage
    {
        const int nx = nx_p[0];
        const int nu = nu_p[0];
        // calculate xIb
        ROWEX(nx - rankI, -alpha_primal, LlIt_p, nx - rankI, 0, delta_ux_p, nu + rankI);
        // assume TRSV_LTN allows aliasing, this is the case in normal BLAS
        TRSV_LTN(nx - rankI, LlIt_p, 0, 0, delta_ux_p, nu + rankI, delta_ux_p, nu + rankI);
        // calculate xIa
        ROWEX(rankI, alpha_primal, GgIt_tilde_p, nx - rankI, 0, delta_ux_p, nu);
        // assume aliasing is possible for last two elements
        GEMV_T(nx - rankI, rankI, 1.0, GgIt_tilde_p, 0, 0, delta_ux_p, nu + rankI, 1.0, delta_ux_p, nu, delta_ux_p, nu);
        //// lag
        // ROWEX(rankI, -alpha_primal, Ppt_p, nx, 0, lam_p, 0);
        // assume aliasing is possible for last two elements
        // GEMV_T(nx, rankI, -1.0, Ppt_p, 0, 0, delta_ux_p, nu, 1.0, lam_p, 0, lam_p, 0);
        // U^-T
        // TRSV_LNN(rankI, HhIt_p, 0, 0, lam_p, 0, lam_p, 0);
        // L^-T
        // TRSV_UNU(rankI, rankI, HhIt_p, 0, 0, lam_p, 0, lam_p, 0);
        // (PlI_p)->PtV(rankI, lam_p, 0);
        (PrI_p)->PtV(rankI, delta_ux_p, nu);
        AXPY(nx, 1.0, delta_ux_p, nu, curr_ux_p, nu, next_ux_p, nu);
    }
    int *offs_ux = (int *)OCP_->aux.ux_offs.data();
    int *offs_g = (int *)OCP_->aux.g_offs.data();
    int *offs_dyn_eq = (int *)OCP_->aux.dyn_eq_offs.data();
    // other stages
    // for (int k = 0; k < K - 1; k++)
    // int dyn_eqs_ofs = offs_g[K - 1] + ng_p[K - 1]; // this value is incremented at end of recursion
    for (int k = 0; k < K - 1; k++)
    {
        const int nx = nx_p[k];
        const int nu = nu_p[k];
        const int nxp1 = nx_p[k + 1];
        const int nup1 = nu_p[k + 1];
        const int offsp1 = offs_ux[k + 1];
        const int offs = offs_ux[k];
        const int rho_k = rho_p[k];
        const int numrho_k = nu - rho_k;
        const int offs_g_k = offs_g[k];
        const int offs_dyn_eq_k = offs_dyn_eq[k];
        const int offs_g_kp1 = offs_g[k + 1];
        const int gammamrho_k = gamma_p[k] - rho_p[k];
        const int gamma_k = gamma_p[k];
        const int gammamrho_kp1 = gamma_p[k + 1] - rho_p[k + 1];
        if (numrho_k > 0)
        {
            /// calculate ukb_tilde
            // -Lkxk - lk
            ROWEX(numrho_k, -alpha_primal, Llt_p + k, numrho_k + nx, 0, delta_ux_p, offs + rho_k);
            // assume aliasing of last two eliments is allowed
            GEMV_T(nx, numrho_k, -1.0, Llt_p + k, numrho_k, 0, delta_ux_p, offs + nu, 1.0, delta_ux_p, offs + rho_k, delta_ux_p, offs + rho_k);
            TRSV_LTN(numrho_k, Llt_p + k, 0, 0, delta_ux_p, offs + rho_k, delta_ux_p, offs + rho_k);
        }
        /// calcualate uka_tilde
        if (rho_k > 0)
        {
            ROWEX(rho_k, alpha_primal, Ggt_tilde_p + k, numrho_k + nx, 0, delta_ux_p, offs);
            GEMV_T(nx + numrho_k, rho_k, 1.0, Ggt_tilde_p + k, 0, 0, delta_ux_p, offs + rho_k, 1.0, delta_ux_p, offs, delta_ux_p, offs);
            // calculate lamda_tilde_k
            // copy vk to right location
            // we implemented a version of vector copy that starts with copy of last element, to avoid aliasing error
            // VECCPR(gammamrho_k, lam_p, offs_g_k, lam_p, offs_g_k + rho_k);
            // ROWEX(rho_k, -alpha_primal, RSQrqt_tilde_p + k, nu + nx, 0, lam_p, offs_g_k);
            // assume aliasing of last two eliments is allowed
            // GEMV_T(nu + nx, rho_k, -1.0, RSQrqt_tilde_p + k, 0, 0, delta_ux_p, offs, 1.0, lam_p, offs_g_k, lam_p, offs_g_k);
            // nu-rank_k+nx,0
            // needless copy because feature not implemented yet in trsv_lnn
            GECP(rho_k, gamma_k, Ggt_tilde_p + k, nu - rho_k + nx + 1, 0, AL_p, 0, 0);
            // U^-T
            // TRSV_LNN(rho_k, AL_p, 0, 0, lam_p, offs_g_k, lam_p, offs_g_k);
            // L^-T
            // TRSV_UNU(rho_k, gamma_k, AL_p, 0, 0, lam_p, offs_g_k, lam_p, offs_g_k);
            // (Pl_p + k)->PtV(rho_k, lam_p, offs_g_k);
            (Pr_p + k)->PtV(rho_k, delta_ux_p, offs);
        }
        // calculate xkp1
        //// multiple shooting
        // ROWEX(nxp1, alpha_primal, BAbt_p + k, nu + nx, 0, delta_ux_p, offsp1 + nup1);
        // GEMV_T(nu + nx, nxp1, 1.0, BAbt_p + k, 0, 0, delta_ux_p, offs, 1.0, delta_ux_p, offsp1 + nup1, delta_ux_p, offsp1 + nup1);
        // AXPY(nxp1, 1.0, delta_ux_p, offsp1+nup1, curr_ux_p, offsp1+nup1, next_ux_p, offsp1+nup1);
        // AXPY(nu, 1.0, delta_ux_p, offs, curr_ux_p, offs, next_ux_p, offs);



        //// DDP
        AXPY(nu, 1.0, delta_ux_p, offs, curr_ux_p, offs, next_ux_p, offs);
        FatropVecBF xp1vec = FatropVecBF(nxp1, offsp1+nup1, next_ux_p);
        ocpinterface_ -> EvalDynamics(OCP_, k, FatropVecBF(nu, offs, next_ux_p), FatropVecBF(nx, offs+nu, next_ux_p), xp1vec);
        ROWEX(nxp1,1.0-alpha_primal, BAbt_p + k, nu + nx, 0, delta_ux_p, offsp1 + nup1);
        AXPY(nxp1, -1.0, delta_ux_p, offsp1+nup1, next_ux_p, offsp1+nup1, next_ux_p, offsp1+nup1);
        AXPY(nxp1, -1.0, curr_ux_p, offsp1+nup1, next_ux_p, offsp1+nup1, delta_ux_p, offsp1+nup1);
        


        // calculate lam_dyn xp1
        // ROWEX(nxp1, alpha_primal, Ppt_p + (k + 1), nxp1, 0, lam_p, offs_dyn_eq_k);
        // GEMV_T(nxp1, nxp1, 1.0, Ppt_p + (k + 1), 0, 0, delta_ux_p, offsp1 + nup1, 1.0, lam_p, offs_dyn_eq_k, lam_p, offs_dyn_eq_k);
        // GEMV_T(gammamrho_kp1, nxp1, 1.0, Hh_p + (k + 1), 0, 0, lam_p, offs_g_kp1, 1.0, lam_p, offs_dyn_eq_k, lam_p, offs_dyn_eq_k);
    }
    for (int k = 0; k < K; k++)
    {
        const int nx = nx_p[k];
        const int nu = nu_p[k];
        const int ng_ineq = ng_ineq_p[k];
        const int offs = offs_ux[k];
        const int offs_g_ineq_k = offs_g_ineq_p[k];
        const int offs_ineq_k = offs_ineq_p[k];
        if (ng_ineq > 0)
        {
            // calculate delta_s
            ROWEX(ng_ineq, alpha_primal, Ggt_ineq_p + k, nu + nx, 0, delta_s_p, offs_ineq_k);
            // GEMV_T(nu + nx, ng_ineq, 1.0, Ggt_ineq_p + k, 0, 0, ux_p, offs, 1.0, delta_s_p, offs_ineq_k, delta_s_p, offs_ineq_k);
            GEMV_T(nu + nx, ng_ineq, 1.0, Ggt_ineq_p + k, 0, 0, delta_ux_p, offs, 1.0, delta_s_p, offs_ineq_k, delta_s_p, offs_ineq_k);
            // calculate lamineq
            for (int i = 0; i < ng_ineq; i++)
            {
                double scaling_factor_L = 0.0;
                double scaling_factor_U = 0.0;
                double zLi = VECEL(zL_p, offs_ineq_k + i);
                double zUi = VECEL(zU_p, offs_ineq_k + i);
                double si = VECEL(s_p, offs_ineq_k + i);
                double loweri = VECEL(lower_p, offs_ineq_k + i);
                double upperi = VECEL(upper_p, offs_ineq_k + i);
                double grad_barrier_L = 0.0;
                double grad_barrier_U = 0.0;
                bool lower_bounded = !isinf(loweri);
                bool upper_bounded = !isinf(upperi);
                double ds = VECEL(delta_s_p, offs_ineq_k + i);
                double lamIi = inertia_correction * ds;
                if (lower_bounded)
                {
                    double dist = si - loweri;
                    double dist_m1 = 1.0 / dist;
                    scaling_factor_L = zLi * dist_m1;
                    grad_barrier_L = -mu * dist_m1;
                    double z = VECEL(zL_p, offs_ineq_k + i);
                    double dz = -grad_barrier_L - z - scaling_factor_L * ds;
                    // VECEL(delta_zL_p, offs_ineq_k + i) = dz;
                    lamIi += -z - dz;
                }
                if (upper_bounded)
                {
                    double dist = upperi - si;
                    double dist_m1 = 1.0 / dist;
                    scaling_factor_U = zUi * dist_m1;
                    grad_barrier_U = mu * dist_m1;
                    double z = VECEL(zU_p, offs_ineq_k + i);
                    double dz = grad_barrier_U - z + scaling_factor_U * ds;
                    // VECEL(delta_zU_p, offs_ineq_k + i) = dz;
                    lamIi += +z + dz;
                    // VECEL(delta_zU_p, offs_ineq_k + i) = dz;
                    // VECEL(delta_zU_p, offs_ineq_k + i) = grad_barrier_U - VECEL(zU_p, offs_ineq_k + i) + scaling_factor_U * VECEL(delta_s_p, offs_ineq_k + i);
                }
                // double grad_barrier = grad_barrier_L + grad_barrier_U;
                if (!(lower_bounded && upper_bounded))
                {
                    lamIi += lower_bounded ? kappa_d * mu : -kappa_d * mu;
                }
                // VECEL(lam_p, offs_g_ineq_k + i) = lamIi - VECEL(lam_curr_p, offs_g_ineq_k + i);
                // VECEL(lam_p, offs_g_ineq_k + i) = grad_barrier + (inertia_correction + scaling_factor_L + scaling_factor_U) * VECEL(delta_s_p, offs_ineq_k + i);
            }
        }
    }
    double alpha_max_pr = 0.0;
    double alpha_max_du = 0.0;
    fatropdata_->AlphaMax(alpha_max_pr, alpha_max_du, MAX(1 - mu, 0.99));
    // cout << "mu " << mu << endl;
    // cout << "alpha prim " << alpha_primal << endl;
    // cout << "alpha dual " << alpha_dual << endl;
    // cout << "alpha prim " << alpha_max_pr << endl;
    // cout << "alpha dual " << alpha_max_du << endl;
    // axpy(alpha_max_pr,fatropdata_-> delta_x,fatropdata_-> x_curr,fatropdata_-> x_next);
    axpy(alpha_max_pr,fatropdata_-> delta_s, fatropdata_->s_curr,fatropdata_-> s_next);
    axpy(alpha_dual,fatropdata_-> delta_zL,fatropdata_-> zL_curr,fatropdata_-> zL_next);
    axpy(alpha_dual,fatropdata_-> delta_zU,fatropdata_-> zU_curr,fatropdata_-> zU_next);
    axpy(alpha_primal,fatropdata_-> lam_calc,fatropdata_-> lam_curr,fatropdata_-> lam_next);
    // axpby(alpha_primal, lam_calc, 1.0 - alpha_primal, lam_curr, lam_next);
    // reset evaluation flags
    fatropdata_->cache_next = FatropData::EvalCache();
    return 0;
}
