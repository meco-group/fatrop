#ifndef OCPLSRICCATIINCLUDED
#define OCPLSRICCATIINCLUDED
#include "OCPKKT.hpp"
#include "OCPLinearSolver.hpp"
#include <cmath>
namespace fatrop
{
    bool check_reg(const int m, MAT *sA, const int ai, const int aj)
    {
        for (int i = 0; i < m; i++)
        {
            if (MATEL(sA, ai + i, aj + i) < 1e-12)
                return false;
        }
        return true;
    }
    class OCPLSRiccati : public OCPLinearSolver
    {
    public:
        OCPLSRiccati(const OCPDims &dims) : Ppt(dims.nx + 1, dims.nx, dims.K),
                                            Hh(dims.nx, dims.nx + 1, dims.K), // the number of eqs can never exceed nx
                                            AL(vector<int>(1, max(dims.nu + dims.nx + 1)), vector<int>(1, max(dims.nx)), 1),
                                            RSQrqt_tilde(dims.nu + dims.nx + 1, dims.nx + dims.nu, dims.K), // TODO, only save first rho rows (can never exceed nu)
                                            Ggt_stripe(vector<int>(1, max(dims.nu + dims.nx + 1)), vector<int>(1, max(dims.nx + dims.nu)), 1),
                                            Ggt_tilde(dims.nu + dims.nx + 1, dims.nx + dims.nu, dims.K), // TODO, only save first rho rows (can never exceed nu)
                                            GgLt(vector<int>(1, max(dims.nu + dims.nx + 1)), vector<int>(1, max(dims.nu + dims.nx)), 1),
                                            RSQrqt_hat(vector<int>(1, max(dims.nu + dims.nx + 1)), vector<int>(1, max(dims.nx + dims.nu)), 1),
                                            Llt(dims.nu + dims.nx + 1, dims.nu, dims.K),
                                            Llt_shift(vector<int>(1, max(dims.nu + dims.nx + 1)), vector<int>(1, max(dims.nu)), 1),
                                            GgIt_tilde(vector<int>(1, dims.nx.at(0) + 1), vector<int>(1, dims.nx.at(0)), 1),
                                            GgLIt(vector<int>(1, dims.nx.at(0) + 1), vector<int>(1, dims.nx.at(0)), 1),
                                            HhIt(vector<int>(1, dims.nx.at(0) + 1), vector<int>(1, dims.nx.at(0)), 1),
                                            PpIt_hat(vector<int>(1, dims.nx.at(0) + 1), vector<int>(1, dims.nx.at(0)), 1),
                                            LlIt(vector<int>(1, dims.nx.at(0) + 1), vector<int>(1, dims.nx.at(0)), 1),
                                            Ggt_ineq_temp(vector<int>(1, max(dims.nu + dims.nx) + 1), vector<int>(1, max(dims.ng_ineq)), 1),
                                            Pl(max(dims.nu), dims.K), // number of equations can never exceed nx
                                            Pr(max(dims.nu), dims.K),
                                            PlI(dims.nx.at(0), 1),
                                            PrI(dims.nx.at(0), 1),
                                            gamma(vector<int>(dims.K, 0)),
                                            rho(vector<int>(dims.K, 0)){};
        // solve a KKT system
        int computeSD(
            OCPKKTMemory *OCP,
            const double inertia_correction_w,
            const double inertia_correction_c,
            const double mu,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &s,
            const FatropVecBF &zL,
            const FatropVecBF &zU,
            const FatropVecBF &delta_zL,
            const FatropVecBF &delta_zU,
            const FatropVecBF &lower,
            const FatropVecBF &upper,
            const FatropVecBF &delta_s) override
        {
            if (inertia_correction_c == 0.0)
            {
                return computeSDnor(OCP, inertia_correction_w, mu, ux, lam, s, zL, zU, delta_zL, delta_zU, lower, upper, delta_s);
            }
            else
            {
                return computeSDDeg(OCP, inertia_correction_w, inertia_correction_c, ux, lam, s, zL, zU, lower, upper, delta_s);
            }
        }
        // solve a KKT system
        int computeSDDeg(
            OCPKKTMemory *OCP,
            const double inertia_correction_w,
            const double inertia_correction_c,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &s,
            const FatropVecBF &zL,
            const FatropVecBF &zU,
            const FatropVecBF &lower,
            const FatropVecBF &upper,
            const FatropVecBF &delta_s)
        {
            // blasfeo_timer timer;
            // blasfeo_tic(&timer);
            // define compiler macros for notational convenience
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP->name)
#define AUXMACRO(type, name, suffix) type name##suffix = ((type)OCP->aux.name)
#define SOLVERMACRO(type, name, suffix) type name##suffix = ((type)name)
            int K = OCP->K;
            // make variables local for efficiency
            OCPMACRO(MAT *, RSQrqt, _p);
            OCPMACRO(MAT *, BAbt, _p);
            OCPMACRO(MAT *, Ggt, _p);
            SOLVERMACRO(MAT *, Ppt, _p);
            SOLVERMACRO(MAT *, AL, _p);
            SOLVERMACRO(MAT *, RSQrqt_tilde, _p);
            SOLVERMACRO(MAT *, Ggt_tilde, _p);
            SOLVERMACRO(MAT *, Llt, _p);
            SOLVERMACRO(MAT *, Llt_shift, _p);
            SOLVERMACRO(MAT *, LlIt, _p);
            SOLVERMACRO(VEC *, ux, _p);
            SOLVERMACRO(VEC *, lam, _p);
            OCPMACRO(int *, nu, _p);
            OCPMACRO(int *, nx, _p);
            OCPMACRO(int *, ng, _p);
            MAT *RSQrq_hat_curr_p;
            double delta_cmin1 = 1 / inertia_correction_c;

            /////////////// recursion ///////////////

            // last stage
            {
                const int nx = nx_p[K - 1];
                const int nu = nu_p[K - 1]; // this should be zero but is included here in case of misuse
                const int ng = ng_p[K - 1];
                // Pp_Km1 <- Qq_Km1
                GECP(nx + 1, nx, RSQrqt_p + (K - 1), nu, nu, Ppt_p + K - 1, 0, 0);
                DIARE(nx, inertia_correction_w, Ppt_p + K - 1, 0, 0);
                GECP(nx + 1, nx, Ggt_p + K - 1, nu, 0, Ggt_tilde_p + K - 1, 0, 0); // needless operation because feature not implemented yet
                SYRK_LN_MN(nx + 1, nx, ng, delta_cmin1, Ggt_tilde_p + K - 1, 0, 0, Ggt_tilde_p + K - 1, 0, 0, 1.0, Ppt_p + K - 1, 0, 0, Ppt_p + K - 1, 0, 0);
                TRTR_L(nx, Ppt_p + K - 1, 0, 0, Ppt_p + K - 1, 0, 0);
            }
            for (int k = K - 2; k >= 0; --k)
            {
                const int nu = nu_p[k];
                const int nx = nx_p[k];
                const int nxp1 = nx_p[k + 1];
                const int ng = ng_p[k];
                //////// SUBSDYN
                {
                    // AL <- [BAb]^T_k P_kp1
                    GEMM_NT(nu + nx + 1, nxp1, nxp1, 1.0, BAbt_p + k, 0, 0, Ppt_p + k + 1, 0, 0, 0.0, AL_p, 0, 0, AL_p, 0, 0);
                    // AL[-1,:] <- AL[-1,:] + p_kp1^T
                    GEAD(1, nxp1, 1.0, Ppt_p + (k + 1), nxp1, 0, AL_p, nx + nu, 0);
                    // RSQrqt_stripe <- AL[BA] + RSQrqt
                    SYRK_LN_MN(nu + nx + 1, nu + nx, nxp1, 1.0, AL_p, 0, 0, BAbt_p + k, 0, 0, 1.0, RSQrqt_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
                    // RSQrqt + 1/d_c At@A
                    DIARE(nu + nx, inertia_correction_w, RSQrqt_tilde_p + k, 0, 0);
                    SYRK_LN_MN(nu + nx + 1, nu + nx, ng, delta_cmin1, Ggt_p + k, 0, 0, Ggt_p + k, 0, 0, 1.0, RSQrqt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
                }
                //////// TRANSFORM_AND_SUBSEQ
                {
                    RSQrq_hat_curr_p = RSQrqt_tilde_p + k;
                }
                //////// SCHUR
                {
                    // DLlt_k = [chol(R_hatk); Llk@chol(R_hatk)^-T]
                    POTRF_L_MN(nu + nx + 1, nu, RSQrq_hat_curr_p, 0, 0, Llt_p + k, 0, 0);
                    if (!check_reg(nu, Llt_p + k, 0, 0))
                        return 1;
                    // Pp_k = Qq_hatk - L_k^T @ Ll_k
                    // SYRK_LN_MN(nx+1, nx, nu-rank_k, -1.0,Llt_p+k, nu-rank_k,0, Llt_p+k, nu-rank_k,0, 1.0, RSQrq_hat_curr_p, nu-rank_k, nu-rank_k,Pp+k,0,0); // feature not implmented yet
                    GECP(nx + 1, nu, Llt_p + k, nu, 0, Llt_shift_p, 0, 0); // needless operation because feature not implemented yet
                    SYRK_LN_MN(nx + 1, nx, nu, -1.0, Llt_shift_p, 0, 0, Llt_shift_p, 0, 0, 1.0, RSQrq_hat_curr_p, nu, nu, Ppt_p + k, 0, 0);
                }
                TRTR_L(nx, Ppt_p + k, 0, 0, Ppt_p + k, 0, 0);
            }
            //////// FIRST_STAGE
            {
                const int nx = nx_p[0];
                {
                    POTRF_L_MN(nx + 1, nx, Ppt_p, 0, 0, LlIt_p, 0, 0);
                    if (!check_reg(nx, LlIt_p, 0, 0))
                        return 2;
                }
            }
            ////// FORWARD_SUBSTITUTION:
            // first stage
            {
                const int nx = nx_p[0];
                const int nu = nu_p[0];
                // calculate xIb
                ROWEX(nx, -1.0, LlIt_p, nx, 0, ux_p, nu);
                // assume TRSV_LTN allows aliasing, this is the case in normal BLAS
                TRSV_LTN(nx, LlIt_p, 0, 0, ux_p, nu, ux_p, nu);
            }
            int *offs_ux = (int *)OCP->aux.ux_offs.data();
            int *offs_g = (int *)OCP->aux.g_offs.data();
            int *offs_dyn_eq = (int *)OCP->aux.dyn_eq_offs.data();
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
                const int offs_dyn_eq_k = offs_dyn_eq[k];
                /// calculate ukb_tilde
                // -Lkxk - lk
                ROWEX(nu, -1.0, Llt_p + k, nu + nx, 0, ux_p, offs);
                // assume aliasing of last two eliments is allowed
                GEMV_T(nx, nu, -1.0, Llt_p + k, nu, 0, ux_p, offs + nu, 1.0, ux_p, offs, ux_p, offs);
                TRSV_LTN(nu, Llt_p + k, 0, 0, ux_p, offs, ux_p, offs);
                // calculate xkp1
                ROWEX(nxp1, 1.0, BAbt_p + k, nu + nx, 0, ux_p, offsp1 + nup1);
                GEMV_T(nu + nx, nxp1, 1.0, BAbt_p + k, 0, 0, ux_p, offs, 1.0, ux_p, offsp1 + nup1, ux_p, offsp1 + nup1);
                // calculate lam_dyn xp1
                ROWEX(nxp1, 1.0, Ppt_p + (k + 1), nxp1, 0, lam_p, offs_dyn_eq_k);
                GEMV_T(nxp1, nxp1, 1.0, Ppt_p + (k + 1), 0, 0, ux_p, offsp1 + nup1, 1.0, lam_p, offs_dyn_eq_k, lam_p, offs_dyn_eq_k);
                // // calculate lam_eq xk
                // ROWEX(ng, -1.0, Ggt_p + k, 0, 0, lam_p, offs_g_k);
                // GEMV_T(nu + nx, ng, delta_cmin1, Ggt_p + k, 0, 0, ux_p, offs, 1.0, lam_p, offs_g_k, lam_p, offs_g_k);
            }
            // calculate lam_eq xk
            for (int k = 0; k < K; k++)
            {
                const int nx = nx_p[k];
                const int nu = nu_p[k];
                const int ng = ng_p[k];
                const int offs = offs_ux[k];
                const int offs_g_k = offs_g[k];
                ROWEX(ng, delta_cmin1, Ggt_p + k, nu + nx, 0, lam_p, offs_g_k);
                GEMV_T(nu + nx, ng, delta_cmin1, Ggt_p + k, 0, 0, ux_p, offs, 1.0, lam_p, offs_g_k, lam_p, offs_g_k);
            }
            // double el = blasfeo_toc(&timer);
            // cout << "el time " << el << endl;
            return 0;
        }
        // solve a KKT system
        int
        SolveInitialization(
            OCPKKTMemory *OCP,
            const FatropVecBF &lam,
            const FatropVecBF &ux_dummy,
            const FatropVecBF &s_dummy,
            const FatropVecBF &zL,
            const FatropVecBF &zU,
            const FatropVecBF &lower,
            const FatropVecBF &upper) override
        {
            // blasfeo_timer timer;
            // blasfeo_tic(&timer);
            // define compiler macros for notational convenience
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP->name)
#define AUXMACRO(type, name, suffix) type name##suffix = ((type)OCP->aux.name)
#define SOLVERMACRO(type, name, suffix) type name##suffix = ((type)name)
            int K = OCP->K;
            // make variables local for efficiency
            OCPMACRO(MAT *, RSQrqt, _p);
            OCPMACRO(MAT *, BAbt, _p);
            OCPMACRO(MAT *, Ggt, _p);
            SOLVERMACRO(MAT *, Ppt, _p);
            SOLVERMACRO(MAT *, Hh, _p);
            SOLVERMACRO(MAT *, AL, _p);
            SOLVERMACRO(MAT *, RSQrqt_tilde, _p);
            SOLVERMACRO(MAT *, Ggt_stripe, _p);
            SOLVERMACRO(MAT *, Ggt_tilde, _p);
            SOLVERMACRO(PMAT *, Pl, _p);
            SOLVERMACRO(PMAT *, Pr, _p);
            SOLVERMACRO(MAT *, GgLt, _p);
            SOLVERMACRO(MAT *, RSQrqt_hat, _p);
            SOLVERMACRO(MAT *, Llt, _p);
            SOLVERMACRO(MAT *, Llt_shift, _p);
            SOLVERMACRO(MAT *, GgIt_tilde, _p);
            SOLVERMACRO(MAT *, GgLIt, _p);
            SOLVERMACRO(MAT *, HhIt, _p);
            SOLVERMACRO(MAT *, PpIt_hat, _p);
            SOLVERMACRO(MAT *, LlIt, _p);
            SOLVERMACRO(PMAT *, PlI, _p);
            SOLVERMACRO(PMAT *, PrI, _p);
            SOLVERMACRO(VEC *, ux_dummy, _p);
            SOLVERMACRO(VEC *, lam, _p);
            OCPMACRO(int *, nu, _p);
            OCPMACRO(int *, nx, _p);
            OCPMACRO(int *, ng, _p);
            SOLVERMACRO(int *, gamma, _p);
            SOLVERMACRO(int *, rho, _p);
            MAT *RSQrq_hat_curr_p;
            int rank_k;

            /////////////// recursion ///////////////

            // last stage
            {
                const int nx = nx_p[K - 1];
                const int nu = nu_p[K - 1]; // this should be zero but is included here in case of misuse
                const int ng = ng_p[K - 1];
                // Pp_Km1 <- Qq_Km1
                GECP(nx + 1, nx, RSQrqt_p + (K - 1), nu, nu, Ppt_p + K - 1, 0, 0);
                // DIARE(nx, inertia_correction, Ppt_p + K - 1, 0, 0);
                // Hh_Km1 <- Gg_Km1
                GETR(nx + 1, ng, Ggt_p + (K - 1), nu, 0, Hh_p + (K - 1), 0, 0);
                gamma_p[K - 1] = ng;
                rho_p[K - 1] = 0;
            }
            for (int k = K - 2; k >= 0; --k)
            {
                const int nu = nu_p[k];
                const int nx = nx_p[k];
                const int nxp1 = nx_p[k + 1];
                const int ng = ng_p[k];
                // calculate the size of H_{k+1} matrix
                const int Hp1_size = gamma_p[k + 1] - rho_p[k + 1];
                if (Hp1_size > nu + nx)
                    return -1;
                // gamma_k <- number of eqs represented by Ggt_stripe
                const int gamma_k = Hp1_size + ng;
                //////// SUBSDYN
                {
                    // AL <- [BAb]^T_k P_kp1
                    GEMM_NT(nu + nx + 1, nxp1, nxp1, 1.0, BAbt_p + k, 0, 0, Ppt_p + k + 1, 0, 0, 0.0, AL_p, 0, 0, AL_p, 0, 0);
                    // AL[-1,:] <- AL[-1,:] + p_kp1^T
                    GEAD(1, nxp1, 1.0, Ppt_p + (k + 1), nxp1, 0, AL_p, nx + nu, 0);
                    // RSQrqt_stripe <- AL[BA] + RSQrqt
                    SYRK_LN_MN(nu + nx + 1, nu + nx, nxp1, 1.0, AL_p, 0, 0, BAbt_p + k, 0, 0, 1.0, RSQrqt_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
                    // DIARE(nu + nx, inertia_correction, RSQrqt_tilde_p + k, 0, 0);
                    gamma_p[k] = gamma_k;
                    // if ng[k]>0
                    if (gamma_k > 0)
                    {
                        // if Gk nonempty
                        if (ng > 0)
                        {
                            // Ggt_stripe  <- Ggt_k
                            GECP(nu + nx + 1, ng, Ggt_p + k, 0, 0, Ggt_stripe_p, 0, 0);
                        }
                        // if Hkp1 nonempty
                        if (Hp1_size > 0)
                        {
                            // Ggt_stripe <- [Ggt_k [BAb_k^T]H_kp1]
                            GEMM_NT(nu + nx + 1, Hp1_size, nxp1, 1.0, BAbt_p + k, 0, 0, Hh_p + (k + 1), 0, 0, 0.0, Ggt_stripe_p, 0, ng, Ggt_stripe_p, 0, ng);
                            // Ggt_stripe[-1,ng:] <- Ggt_stripe[-1,ng:] + h_kp1^T
                            GEADTR(1, Hp1_size, 1.0, Hh_p + (k + 1), 0, nxp1, Ggt_stripe_p, nu + nx, ng);
                        }
                    }
                    else
                    {
                        rho_p[k] = 0;
                        rank_k = 0;
                        RSQrq_hat_curr_p = RSQrqt_tilde_p + k;
                    }
                }
                //////// TRANSFORM_AND_SUBSEQ
                {
                    // symmetric transformation, done a little different than in paper, in order to fuse LA operations
                    // LU_FACT_TRANSPOSE(Ggtstripe[:gamma_k, nu+nx+1], nu max)
                    LU_FACT_transposed(gamma_k, nu + nx + 1, nu, rank_k, Ggt_stripe_p, Pl_p + k, Pr_p + k);
                    rho_p[k] = rank_k;
                    if (gamma_k - rank_k > 0)
                    {
                        // transfer eq's to next stage
                        GETR(nx + 1, gamma_k - rank_k, Ggt_stripe_p, nu, rank_k, Hh_p + k, 0, 0);
                    }
                    if (rank_k > 0)
                    {
                        // Ggt_tilde_k <- Ggt_stripe[rho_k:nu+nx+1, :rho] L-T (note that this is slightly different from the implementation)
                        TRSM_RLNN(nu - rank_k + nx + 1, rank_k, -1.0, Ggt_stripe_p, 0, 0, Ggt_stripe_p, rank_k, 0, Ggt_tilde_p + k, 0, 0);
                        // the following command copies the top block matrix (LU) to the bottom because it it needed later
                        GECP(rank_k, gamma_k, Ggt_stripe_p, 0, 0, Ggt_tilde_p + k, nu - rank_k + nx + 1, 0);
                        // permutations
                        TRTR_L(nu + nx, RSQrqt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0); // copy lower part of RSQ to upper part
                        (Pr_p + k)->PM(rank_k, RSQrqt_tilde_p + k);                          // TODO make use of symmetry
                        (Pr_p + k)->MPt(rank_k, RSQrqt_tilde_p + k);
                        // GL <- Ggt_tilde_k @ RSQ[:rho,:nu+nx] + RSQrqt[rho:nu+nx+1, rho:] (with RSQ[:rho,:nu+nx] = RSQrqt[:nu+nx,:rho]^T)
                        // GEMM_NT(nu - rank_k + nx + 1, nu + nx, rank_k, 1.0, Ggt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0, 1.0, RSQrqt_tilde_p + k, rank_k, 0, GgLt_p, 0, 0);
                        // split up because valgrind was giving invalid read errors when C matrix has nonzero row offset
                        GECP(nu - rank_k + nx + 1, nu + nx, RSQrqt_tilde_p + k, rank_k, 0, GgLt_p, 0, 0);
                        GEMM_NT(nu - rank_k + nx + 1, nu + nx, rank_k, 1.0, Ggt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0, 1.0, GgLt_p, 0, 0, GgLt_p, 0, 0);
                        // RSQrqt_hat = GgLt[nu-rank_k + nx +1, :rank_k] * G[:rank_k, :nu+nx] + GgLt[rank_k:, :]  (with G[:rank_k,:nu+nx] = Gt[:nu+nx,:rank_k]^T)
                        SYRK_LN_MN(nu - rank_k + nx + 1, nu + nx - rank_k, rank_k, 1.0, GgLt_p, 0, 0, Ggt_tilde_p + k, 0, 0, 1.0, GgLt_p, 0, rank_k, RSQrqt_hat_p, 0, 0);
                        // GEMM_NT(nu - rank_k + nx + 1, nu + nx - rank_k, rank_k, 1.0, GgLt_p, 0, 0, Ggt_tilde_p + k, 0, 0, 1.0, GgLt_p, 0, rank_k, RSQrqt_hat_p, 0, 0);
                        RSQrq_hat_curr_p = RSQrqt_hat_p;
                    }
                    else
                    {
                        RSQrq_hat_curr_p = RSQrqt_tilde_p + k;
                    }
                }
                //////// SCHUR
                {
                    if (nu - rank_k > 0)
                    {
                        // DLlt_k = [chol(R_hatk); Llk@chol(R_hatk)^-T]
                        POTRF_L_MN(nu - rank_k + nx + 1, nu - rank_k, RSQrq_hat_curr_p, 0, 0, Llt_p + k, 0, 0);
                        if (!check_reg(nu - rank_k, Llt_p + k, 0, 0))
                            return 1;
                        // Pp_k = Qq_hatk - L_k^T @ Ll_k
                        // SYRK_LN_MN(nx+1, nx, nu-rank_k, -1.0,Llt_p+k, nu-rank_k,0, Llt_p+k, nu-rank_k,0, 1.0, RSQrq_hat_curr_p, nu-rank_k, nu-rank_k,Pp+k,0,0); // feature not implmented yet
                        GECP(nx + 1, nu - rank_k, Llt_p + k, nu - rank_k, 0, Llt_shift_p, 0, 0); // needless operation because feature not implemented yet
                        // SYRK_LN_MN(nx + 1, nx, nu - rank_k, -1.0, Llt_shift_p, 0, 0, Llt_shift_p, 0, 0, 1.0, RSQrq_hat_curr_p, nu - rank_k, nu - rank_k, Ppt_p + k, 0, 0);
                        GECP(nx + 1, nx, RSQrq_hat_curr_p, nu - rank_k, nu - rank_k, Ppt_p + k, 0, 0);
                        SYRK_LN_MN(nx + 1, nx, nu - rank_k, -1.0, Llt_shift_p, 0, 0, Llt_shift_p, 0, 0, 1.0, Ppt_p + k, 0, 0, Ppt_p + k, 0, 0);
                    }
                    else
                    {
                        GECP(nx + 1, nx, RSQrq_hat_curr_p, 0, 0, Ppt_p + k, 0, 0);
                    }
                    TRTR_L(nx, Ppt_p + k, 0, 0, Ppt_p + k, 0, 0);
                }
            }
            int rankI = 0;
            //////// FIRST_STAGE
            {
                const int nx = nx_p[0];
                int gamma_I = gamma_p[0] - rho_p[0];
                if (gamma_I > 0)
                {
                    GETR(gamma_I, nx + 1, Hh_p + 0, 0, 0, HhIt_p, 0, 0); // transposition may be avoided
                    // HhIt[0].print();
                    LU_FACT_transposed(gamma_I, nx + 1, nx, rankI, HhIt_p, PlI_p, PrI_p);
                    if (rankI < gamma_I)
                        return -2;
                    // PpIt_tilde <- Ggt[rankI:nx+1, :rankI] L-T (note that this is slightly different from the implementation)
                    TRSM_RLNN(nx - rankI + 1, rankI, -1.0, HhIt_p, 0, 0, HhIt_p, rankI, 0, GgIt_tilde_p, 0, 0);
                    // permutations
                    (PrI_p)->PM(rankI, Ppt_p); // TODO make use of symmetry
                    (PrI_p)->MPt(rankI, Ppt_p);
                    // // GL <- GgIt_tilde @ Pp[:rankI,:nx] + Ppt[rankI:nx+1, rankI:] (with Pp[:rankI,:nx] = Ppt[:nx,:rankI]^T)
                    // GEMM_NT(nx - rankI + 1, nx, rankI, 1.0, GgIt_tilde_p, 0, 0, Ppt_p, 0, 0, 1.0, Ppt_p, rankI, 0, GgLIt_p, 0, 0);
                    // split up because valgrind was giving invalid read errors when C matrix has nonzero row offset
                    GECP(nx - rankI + 1, nx, Ppt_p, rankI, 0, GgLIt_p, 0, 0);
                    GEMM_NT(nx - rankI + 1, nx, rankI, 1.0, GgIt_tilde_p, 0, 0, Ppt_p, 0, 0, 1.0, GgLIt_p, 0, 0, GgLIt_p, 0, 0);
                    // // RSQrqt_hat = GgLt[nu-rank_k + nx +1, :rank_k] * G[:rank_k, :nu+nx] + GgLt[rank_k:, :]  (with G[:rank_k,:nu+nx] = Gt[:nu+nx,:rank_k]^T)
                    SYRK_LN_MN(nx - rankI + 1, nx - rankI, rankI, 1.0, GgLIt_p, 0, 0, GgIt_tilde_p, 0, 0, 1.0, GgLIt_p, 0, rankI, PpIt_hat_p, 0, 0);
                    // TODO skipped if nx-rankI = 0
                    POTRF_L_MN(nx - rankI + 1, nx - rankI, PpIt_hat_p, 0, 0, LlIt_p, 0, 0);
                    if (!check_reg(nx - rankI, LlIt_p, 0, 0))
                        return 2;
                }
                else
                {
                    rankI = 0;
                    POTRF_L_MN(nx + 1, nx, Ppt_p, 0, 0, LlIt_p, 0, 0);
                    if (!check_reg(nx, LlIt_p, 0, 0))
                        return 2;
                }
            }
            ////// FORWARD_SUBSTITUTION:
            // first stage
            {
                const int nx = nx_p[0];
                const int nu = nu_p[0];
                // calculate xIb
                ROWEX(nx - rankI, -1.0, LlIt_p, nx - rankI, 0, ux_dummy_p, nu + rankI);
                // assume TRSV_LTN allows aliasing, this is the case in normal BLAS
                TRSV_LTN(nx - rankI, LlIt_p, 0, 0, ux_dummy_p, nu + rankI, ux_dummy_p, nu + rankI);
                // calculate xIa
                ROWEX(rankI, 1.0, GgIt_tilde_p, nx - rankI, 0, ux_dummy_p, nu);
                // assume aliasing is possible for last two elements
                GEMV_T(nx - rankI, rankI, 1.0, GgIt_tilde_p, 0, 0, ux_dummy_p, nu + rankI, 1.0, ux_dummy_p, nu, ux_dummy_p, nu);
                //// lag
                ROWEX(rankI, -1.0, Ppt_p, nx, 0, lam_p, 0);
                // assume aliasing is possible for last two elements
                GEMV_T(nx, rankI, -1.0, Ppt_p, 0, 0, ux_dummy_p, nu, 1.0, lam_p, 0, lam_p, 0);
                // U^-T
                TRSV_LNN(rankI, HhIt_p, 0, 0, lam_p, 0, lam_p, 0);
                // L^-T
                TRSV_UNU(rankI, rankI, HhIt_p, 0, 0, lam_p, 0, lam_p, 0);
                (PlI_p)->PtV(rankI, lam_p, 0);
                (PrI_p)->PtV(rankI, ux_dummy_p, nu);
            }
            int *offs_ux = (int *)OCP->aux.ux_offs.data();
            int *offs_g = (int *)OCP->aux.g_offs.data();
            int *offs_dyn_eq = (int *)OCP->aux.dyn_eq_offs.data();
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
                    ROWEX(numrho_k, -1.0, Llt_p + k, numrho_k + nx, 0, ux_dummy_p, offs + rho_k);
                    // assume aliasing of last two eliments is allowed
                    GEMV_T(nx, numrho_k, -1.0, Llt_p + k, numrho_k, 0, ux_dummy_p, offs + nu, 1.0, ux_dummy_p, offs + rho_k, ux_dummy_p, offs + rho_k);
                    TRSV_LTN(numrho_k, Llt_p + k, 0, 0, ux_dummy_p, offs + rho_k, ux_dummy_p, offs + rho_k);
                }
                /// calcualate uka_tilde
                if (rho_k > 0)
                {
                    ROWEX(rho_k, 1.0, Ggt_tilde_p + k, numrho_k + nx, 0, ux_dummy_p, offs);
                    GEMV_T(nx + numrho_k, rho_k, 1.0, Ggt_tilde_p + k, 0, 0, ux_dummy_p, offs + rho_k, 1.0, ux_dummy_p, offs, ux_dummy_p, offs);
                    // calculate lamda_tilde_k
                    // copy vk to right location
                    // we implemented a version of vector copy that starts with copy of last element, to avoid aliasing error
                    VECCPR(gammamrho_k, lam_p, offs_g_k, lam_p, offs_g_k + rho_k);
                    ROWEX(rho_k, -1.0, RSQrqt_tilde_p + k, nu + nx, 0, lam_p, offs_g_k);
                    // assume aliasing of last two eliments is allowed
                    GEMV_T(nu + nx, rho_k, -1.0, RSQrqt_tilde_p + k, 0, 0, ux_dummy_p, offs, 1.0, lam_p, offs_g_k, lam_p, offs_g_k);
                    // nu-rank_k+nx,0
                    // needless copy because feature not implemented yet in trsv_lnn
                    GECP(rho_k, gamma_k, Ggt_tilde_p + k, nu - rho_k + nx + 1, 0, AL_p, 0, 0);
                    // U^-T
                    TRSV_LNN(rho_k, AL_p, 0, 0, lam_p, offs_g_k, lam_p, offs_g_k);
                    // L^-T
                    TRSV_UNU(rho_k, gamma_k, AL_p, 0, 0, lam_p, offs_g_k, lam_p, offs_g_k);
                    (Pl_p + k)->PtV(rho_k, lam_p, offs_g_k);
                    (Pr_p + k)->PtV(rho_k, ux_dummy_p, offs);
                }
                // calculate xkp1
                ROWEX(nxp1, 1.0, BAbt_p + k, nu + nx, 0, ux_dummy_p, offsp1 + nup1);
                GEMV_T(nu + nx, nxp1, 1.0, BAbt_p + k, 0, 0, ux_dummy_p, offs, 1.0, ux_dummy_p, offsp1 + nup1, ux_dummy_p, offsp1 + nup1);
                // calculate lam_dyn xp1
                ROWEX(nxp1, 1.0, Ppt_p + (k + 1), nxp1, 0, lam_p, offs_dyn_eq_k);
                GEMV_T(nxp1, nxp1, 1.0, Ppt_p + (k + 1), 0, 0, ux_dummy_p, offsp1 + nup1, 1.0, lam_p, offs_dyn_eq_k, lam_p, offs_dyn_eq_k);
                GEMV_T(gammamrho_kp1, nxp1, 1.0, Hh_p + (k + 1), 0, 0, lam_p, offs_g_kp1, 1.0, lam_p, offs_dyn_eq_k, lam_p, offs_dyn_eq_k);
            }
            // double el = blasfeo_toc(&timer);
            // cout << "el time " << el << endl;
            return 0;
        }
        // solve a KKT system
        int
        computeSDnor(
            OCPKKTMemory *OCP,
            const double inertia_correction,
            const double mu,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &s,
            const FatropVecBF &zL,
            const FatropVecBF &zU,
            const FatropVecBF &delta_zL,
            const FatropVecBF &delta_zU,
            const FatropVecBF &lower,
            const FatropVecBF &upper,
            const FatropVecBF &delta_s)
        {
            // blasfeo_timer timer;
            // blasfeo_tic(&timer);
            // define compiler macros for notational convenience
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP->name)
#define AUXMACRO(type, name, suffix) type name##suffix = ((type)OCP->aux.name)
#define SOLVERMACRO(type, name, suffix) type name##suffix = ((type)name)
            int K = OCP->K;
            // make variables local for efficiency
            OCPMACRO(MAT *, RSQrqt, _p);
            OCPMACRO(MAT *, BAbt, _p);
            OCPMACRO(MAT *, Ggt, _p);
            OCPMACRO(MAT *, Ggt_ineq, _p);
            SOLVERMACRO(MAT *, Ppt, _p);
            SOLVERMACRO(MAT *, Hh, _p);
            SOLVERMACRO(MAT *, AL, _p);
            SOLVERMACRO(MAT *, RSQrqt_tilde, _p);
            SOLVERMACRO(MAT *, Ggt_stripe, _p);
            SOLVERMACRO(MAT *, Ggt_tilde, _p);
            SOLVERMACRO(PMAT *, Pl, _p);
            SOLVERMACRO(PMAT *, Pr, _p);
            SOLVERMACRO(MAT *, GgLt, _p);
            SOLVERMACRO(MAT *, RSQrqt_hat, _p);
            SOLVERMACRO(MAT *, Llt, _p);
            SOLVERMACRO(MAT *, Llt_shift, _p);
            SOLVERMACRO(MAT *, GgIt_tilde, _p);
            SOLVERMACRO(MAT *, GgLIt, _p);
            SOLVERMACRO(MAT *, HhIt, _p);
            SOLVERMACRO(MAT *, PpIt_hat, _p);
            SOLVERMACRO(MAT *, LlIt, _p);
            SOLVERMACRO(MAT *, Ggt_ineq_temp, _p);
            SOLVERMACRO(PMAT *, PlI, _p);
            SOLVERMACRO(PMAT *, PrI, _p);
            SOLVERMACRO(VEC *, ux, _p);
            SOLVERMACRO(VEC *, lam, _p);
            SOLVERMACRO(VEC *, s, _p);
            SOLVERMACRO(VEC *, zL, _p);
            SOLVERMACRO(VEC *, zU, _p);
            SOLVERMACRO(VEC *, delta_zL, _p);
            SOLVERMACRO(VEC *, delta_zU, _p);
            SOLVERMACRO(VEC *, lower, _p);
            SOLVERMACRO(VEC *, upper, _p);
            SOLVERMACRO(VEC *, delta_s, _p);
            OCPMACRO(int *, nu, _p);
            OCPMACRO(int *, nx, _p);
            OCPMACRO(int *, ng, _p);
            OCPMACRO(int *, ng_ineq, _p);
            SOLVERMACRO(int *, gamma, _p);
            SOLVERMACRO(int *, rho, _p);
            MAT *RSQrq_hat_curr_p;
            int rank_k;
            int *offs_ineq_p = (int *)OCP->aux.ineq_offs.data();

            /////////////// recursion ///////////////

            // last stage
            {
                const int nx = nx_p[K - 1];
                const int nu = nu_p[K - 1]; // this should be zero but is included here in case of misuse
                const int ng = ng_p[K - 1];
                // Pp_Km1 <- Qq_Km1
                GECP(nx + 1, nx, RSQrqt_p + (K - 1), nu, nu, Ppt_p + K - 1, 0, 0);
                DIARE(nx, inertia_correction, Ppt_p + K - 1, 0, 0);
                // Hh_Km1 <- Gg_Km1
                GETR(nx + 1, ng, Ggt_p + (K - 1), nu, 0, Hh_p + (K - 1), 0, 0);
                gamma_p[K - 1] = ng;
                rho_p[K - 1] = 0;
            }
            for (int k = K - 2; k >= 0; --k)
            {
                const int nu = nu_p[k];
                const int nx = nx_p[k];
                const int nxp1 = nx_p[k + 1];
                const int ng = ng_p[k];
                const int ng_ineq = ng_ineq_p[k];
                // const int offs_g_k = offs_g_p[k];
                const int offs_ineq_k = offs_ineq_p[k];
                // calculate the size of H_{k+1} matrix
                const int Hp1_size = gamma_p[k + 1] - rho_p[k + 1];
                if (Hp1_size > nu + nx)
                    return -1;
                // gamma_k <- number of eqs represented by Ggt_stripe
                const int gamma_k = Hp1_size + ng;
                //////// SUBSDYN
                {
                    // AL <- [BAb]^T_k P_kp1
                    GEMM_NT(nu + nx + 1, nxp1, nxp1, 1.0, BAbt_p + k, 0, 0, Ppt_p + k + 1, 0, 0, 0.0, AL_p, 0, 0, AL_p, 0, 0);
                    // AL[-1,:] <- AL[-1,:] + p_kp1^T
                    GEAD(1, nxp1, 1.0, Ppt_p + (k + 1), nxp1, 0, AL_p, nx + nu, 0);
                    // RSQrqt_stripe <- AL[BA] + RSQrqt
                    SYRK_LN_MN(nu + nx + 1, nu + nx, nxp1, 1.0, AL_p, 0, 0, BAbt_p + k, 0, 0, 1.0, RSQrqt_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
                    //// inequalities
                    if (ng_ineq > 0)
                    {
                        GECP(nu + nx, ng_ineq, Ggt_ineq_p + k, 0, 0, Ggt_ineq_temp_p, 0, 0);
                        for (int i = 0; i < ng_ineq; i++)
                        {
                            double scaling_factor = inertia_correction;
                            double zLi = VECEL(zL_p, offs_ineq_k + i);
                            double zUi = VECEL(zU_p, offs_ineq_k + i);
                            double si = VECEL(s_p, offs_ineq_k + i);
                            double loweri = VECEL(lower_p, offs_ineq_k + i);
                            double upperi = VECEL(upper_p, offs_ineq_k + i);
                            double grad_barrier = 0.0;
                            if (!isinf(loweri))
                            {
                                double dist = si - loweri;
                                double dist_m1 = 1.0 / dist;
                                scaling_factor += zLi * dist_m1;
                                grad_barrier += mu * dist_m1;
                            }
                            if (!isinf(upperi))
                            {
                                double dist = upperi - si;
                                double dist_m1 = 1.0 / dist;
                                scaling_factor += zUi * dist_m1;
                                grad_barrier -= mu * dist_m1;
                            }
                            COLSC(nu + nx + 1, scaling_factor, Ggt_ineq_temp_p, 0, i);
                            MATEL(Ggt_ineq_temp_p, nu + nx, i) += grad_barrier;
                        }
                        // add the penalty
                        SYRK_LN_MN(nu + nx + 1, nu + nx, ng_ineq, 1.0, Ggt_ineq_temp_p, 0, 0, Ggt_ineq_p + k, 0, 0, 1.0, RSQrqt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
                    }
                    DIARE(nu + nx, inertia_correction, RSQrqt_tilde_p + k, 0, 0);
                    gamma_p[k] = gamma_k;
                    // if ng[k]>0
                    if (gamma_k > 0)
                    {
                        // if Gk nonempty
                        if (ng > 0)
                        {
                            // Ggt_stripe  <- Ggt_k
                            GECP(nu + nx + 1, ng, Ggt_p + k, 0, 0, Ggt_stripe_p, 0, 0);
                        }
                        // if Hkp1 nonempty
                        if (Hp1_size > 0)
                        {
                            // Ggt_stripe <- [Ggt_k [BAb_k^T]H_kp1]
                            GEMM_NT(nu + nx + 1, Hp1_size, nxp1, 1.0, BAbt_p + k, 0, 0, Hh_p + (k + 1), 0, 0, 0.0, Ggt_stripe_p, 0, ng, Ggt_stripe_p, 0, ng);
                            // Ggt_stripe[-1,ng:] <- Ggt_stripe[-1,ng:] + h_kp1^T
                            GEADTR(1, Hp1_size, 1.0, Hh_p + (k + 1), 0, nxp1, Ggt_stripe_p, nu + nx, ng);
                        }
                    }
                    else
                    {
                        rho_p[k] = 0;
                        rank_k = 0;
                        RSQrq_hat_curr_p = RSQrqt_tilde_p + k;
                    }
                }
                //////// TRANSFORM_AND_SUBSEQ
                {
                    // symmetric transformation, done a little different than in paper, in order to fuse LA operations
                    // LU_FACT_TRANSPOSE(Ggtstripe[:gamma_k, nu+nx+1], nu max)
                    LU_FACT_transposed(gamma_k, nu + nx + 1, nu, rank_k, Ggt_stripe_p, Pl_p + k, Pr_p + k);
                    rho_p[k] = rank_k;
                    if (gamma_k - rank_k > 0)
                    {
                        // transfer eq's to next stage
                        GETR(nx + 1, gamma_k - rank_k, Ggt_stripe_p, nu, rank_k, Hh_p + k, 0, 0);
                    }
                    if (rank_k > 0)
                    {
                        // Ggt_tilde_k <- Ggt_stripe[rho_k:nu+nx+1, :rho] L-T (note that this is slightly different from the implementation)
                        TRSM_RLNN(nu - rank_k + nx + 1, rank_k, -1.0, Ggt_stripe_p, 0, 0, Ggt_stripe_p, rank_k, 0, Ggt_tilde_p + k, 0, 0);
                        // the following command copies the top block matrix (LU) to the bottom because it it needed later
                        GECP(rank_k, gamma_k, Ggt_stripe_p, 0, 0, Ggt_tilde_p + k, nu - rank_k + nx + 1, 0);
                        // permutations
                        TRTR_L(nu + nx, RSQrqt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0); // copy lower part of RSQ to upper part
                        (Pr_p + k)->PM(rank_k, RSQrqt_tilde_p + k);                          // TODO make use of symmetry
                        (Pr_p + k)->MPt(rank_k, RSQrqt_tilde_p + k);
                        // GL <- Ggt_tilde_k @ RSQ[:rho,:nu+nx] + RSQrqt[rho:nu+nx+1, rho:] (with RSQ[:rho,:nu+nx] = RSQrqt[:nu+nx,:rho]^T)
                        // GEMM_NT(nu - rank_k + nx + 1, nu + nx, rank_k, 1.0, Ggt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0, 1.0, RSQrqt_tilde_p + k, rank_k, 0, GgLt_p, 0, 0);
                        // split up because valgrind was giving invalid read errors when C matrix has nonzero row offset
                        GECP(nu - rank_k + nx + 1, nu + nx, RSQrqt_tilde_p + k, rank_k, 0, GgLt_p, 0, 0);
                        GEMM_NT(nu - rank_k + nx + 1, nu + nx, rank_k, 1.0, Ggt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0, 1.0, GgLt_p, 0, 0, GgLt_p, 0, 0);
                        // RSQrqt_hat = GgLt[nu-rank_k + nx +1, :rank_k] * G[:rank_k, :nu+nx] + GgLt[rank_k:, :]  (with G[:rank_k,:nu+nx] = Gt[:nu+nx,:rank_k]^T)
                        SYRK_LN_MN(nu - rank_k + nx + 1, nu + nx - rank_k, rank_k, 1.0, GgLt_p, 0, 0, Ggt_tilde_p + k, 0, 0, 1.0, GgLt_p, 0, rank_k, RSQrqt_hat_p, 0, 0);
                        // GEMM_NT(nu - rank_k + nx + 1, nu + nx - rank_k, rank_k, 1.0, GgLt_p, 0, 0, Ggt_tilde_p + k, 0, 0, 1.0, GgLt_p, 0, rank_k, RSQrqt_hat_p, 0, 0);
                        RSQrq_hat_curr_p = RSQrqt_hat_p;
                    }
                    else
                    {
                        RSQrq_hat_curr_p = RSQrqt_tilde_p + k;
                    }
                }
                //////// SCHUR
                {
                    if (nu - rank_k > 0)
                    {
                        // DLlt_k = [chol(R_hatk); Llk@chol(R_hatk)^-T]
                        POTRF_L_MN(nu - rank_k + nx + 1, nu - rank_k, RSQrq_hat_curr_p, 0, 0, Llt_p + k, 0, 0);
                        if (!check_reg(nu - rank_k, Llt_p + k, 0, 0))
                            return 1;
                        // Pp_k = Qq_hatk - L_k^T @ Ll_k
                        // SYRK_LN_MN(nx+1, nx, nu-rank_k, -1.0,Llt_p+k, nu-rank_k,0, Llt_p+k, nu-rank_k,0, 1.0, RSQrq_hat_curr_p, nu-rank_k, nu-rank_k,Pp+k,0,0); // feature not implmented yet
                        GECP(nx + 1, nu - rank_k, Llt_p + k, nu - rank_k, 0, Llt_shift_p, 0, 0); // needless operation because feature not implemented yet
                        // SYRK_LN_MN(nx + 1, nx, nu - rank_k, -1.0, Llt_shift_p, 0, 0, Llt_shift_p, 0, 0, 1.0, RSQrq_hat_curr_p, nu - rank_k, nu - rank_k, Ppt_p + k, 0, 0);
                        GECP(nx + 1, nx, RSQrq_hat_curr_p, nu - rank_k, nu - rank_k, Ppt_p + k, 0, 0);
                        SYRK_LN_MN(nx + 1, nx, nu - rank_k, -1.0, Llt_shift_p, 0, 0, Llt_shift_p, 0, 0, 1.0, Ppt_p + k, 0, 0, Ppt_p + k, 0, 0);
                    }
                    else
                    {
                        GECP(nx + 1, nx, RSQrq_hat_curr_p, 0, 0, Ppt_p + k, 0, 0);
                    }
                    TRTR_L(nx, Ppt_p + k, 0, 0, Ppt_p + k, 0, 0);
                }
            }
            int rankI = 0;
            //////// FIRST_STAGE
            {
                const int nx = nx_p[0];
                int gamma_I = gamma_p[0] - rho_p[0];
                if (gamma_I > 0)
                {
                    GETR(gamma_I, nx + 1, Hh_p + 0, 0, 0, HhIt_p, 0, 0); // transposition may be avoided
                    // HhIt[0].print();
                    LU_FACT_transposed(gamma_I, nx + 1, nx, rankI, HhIt_p, PlI_p, PrI_p);
                    if (rankI < gamma_I)
                        return -2;
                    // PpIt_tilde <- Ggt[rankI:nx+1, :rankI] L-T (note that this is slightly different from the implementation)
                    TRSM_RLNN(nx - rankI + 1, rankI, -1.0, HhIt_p, 0, 0, HhIt_p, rankI, 0, GgIt_tilde_p, 0, 0);
                    // permutations
                    (PrI_p)->PM(rankI, Ppt_p); // TODO make use of symmetry
                    (PrI_p)->MPt(rankI, Ppt_p);
                    // // GL <- GgIt_tilde @ Pp[:rankI,:nx] + Ppt[rankI:nx+1, rankI:] (with Pp[:rankI,:nx] = Ppt[:nx,:rankI]^T)
                    // GEMM_NT(nx - rankI + 1, nx, rankI, 1.0, GgIt_tilde_p, 0, 0, Ppt_p, 0, 0, 1.0, Ppt_p, rankI, 0, GgLIt_p, 0, 0);
                    // split up because valgrind was giving invalid read errors when C matrix has nonzero row offset
                    GECP(nx - rankI + 1, nx, Ppt_p, rankI, 0, GgLIt_p, 0, 0);
                    GEMM_NT(nx - rankI + 1, nx, rankI, 1.0, GgIt_tilde_p, 0, 0, Ppt_p, 0, 0, 1.0, GgLIt_p, 0, 0, GgLIt_p, 0, 0);
                    // // RSQrqt_hat = GgLt[nu-rank_k + nx +1, :rank_k] * G[:rank_k, :nu+nx] + GgLt[rank_k:, :]  (with G[:rank_k,:nu+nx] = Gt[:nu+nx,:rank_k]^T)
                    SYRK_LN_MN(nx - rankI + 1, nx - rankI, rankI, 1.0, GgLIt_p, 0, 0, GgIt_tilde_p, 0, 0, 1.0, GgLIt_p, 0, rankI, PpIt_hat_p, 0, 0);
                    // TODO skipped if nx-rankI = 0
                    POTRF_L_MN(nx - rankI + 1, nx - rankI, PpIt_hat_p, 0, 0, LlIt_p, 0, 0);
                    if (!check_reg(nx - rankI, LlIt_p, 0, 0))
                        return 2;
                }
                else
                {
                    rankI = 0;
                    POTRF_L_MN(nx + 1, nx, Ppt_p, 0, 0, LlIt_p, 0, 0);
                    if (!check_reg(nx, LlIt_p, 0, 0))
                        return 2;
                }
            }
            ////// FORWARD_SUBSTITUTION:
            // first stage
            {
                const int nx = nx_p[0];
                const int nu = nu_p[0];
                // calculate xIb
                ROWEX(nx - rankI, -1.0, LlIt_p, nx - rankI, 0, ux_p, nu + rankI);
                // assume TRSV_LTN allows aliasing, this is the case in normal BLAS
                TRSV_LTN(nx - rankI, LlIt_p, 0, 0, ux_p, nu + rankI, ux_p, nu + rankI);
                // calculate xIa
                ROWEX(rankI, 1.0, GgIt_tilde_p, nx - rankI, 0, ux_p, nu);
                // assume aliasing is possible for last two elements
                GEMV_T(nx - rankI, rankI, 1.0, GgIt_tilde_p, 0, 0, ux_p, nu + rankI, 1.0, ux_p, nu, ux_p, nu);
                //// lag
                ROWEX(rankI, -1.0, Ppt_p, nx, 0, lam_p, 0);
                // assume aliasing is possible for last two elements
                GEMV_T(nx, rankI, -1.0, Ppt_p, 0, 0, ux_p, nu, 1.0, lam_p, 0, lam_p, 0);
                // U^-T
                TRSV_LNN(rankI, HhIt_p, 0, 0, lam_p, 0, lam_p, 0);
                // L^-T
                TRSV_UNU(rankI, rankI, HhIt_p, 0, 0, lam_p, 0, lam_p, 0);
                (PlI_p)->PtV(rankI, lam_p, 0);
                (PrI_p)->PtV(rankI, ux_p, nu);
            }
            int *offs_ux = (int *)OCP->aux.ux_offs.data();
            int *offs_g = (int *)OCP->aux.g_offs.data();
            int *offs_dyn_eq = (int *)OCP->aux.dyn_eq_offs.data();
            int *offs_g_ineq_p = (int *)OCP->aux.g_ineq_offs.data();
            // other stages
            // for (int k = 0; k < K - 1; k++)
            // int dyn_eqs_ofs = offs_g[K - 1] + ng_p[K - 1]; // this value is incremented at end of recursion
            for (int k = 0; k < K - 1; k++)
            {
                const int nx = nx_p[k];
                const int nu = nu_p[k];
                const int nxp1 = nx_p[k + 1];
                const int nup1 = nu_p[k + 1];
                const int ng_ineq = ng_ineq_p[k];
                const int offsp1 = offs_ux[k + 1];
                const int offs = offs_ux[k];
                const int offs_g_ineq_k = offs_g_ineq_p[k];
                const int offs_ineq_k = offs_ineq_p[k];
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
                    ROWEX(numrho_k, -1.0, Llt_p + k, numrho_k + nx, 0, ux_p, offs + rho_k);
                    // assume aliasing of last two eliments is allowed
                    GEMV_T(nx, numrho_k, -1.0, Llt_p + k, numrho_k, 0, ux_p, offs + nu, 1.0, ux_p, offs + rho_k, ux_p, offs + rho_k);
                    TRSV_LTN(numrho_k, Llt_p + k, 0, 0, ux_p, offs + rho_k, ux_p, offs + rho_k);
                }
                /// calcualate uka_tilde
                if (rho_k > 0)
                {
                    ROWEX(rho_k, 1.0, Ggt_tilde_p + k, numrho_k + nx, 0, ux_p, offs);
                    GEMV_T(nx + numrho_k, rho_k, 1.0, Ggt_tilde_p + k, 0, 0, ux_p, offs + rho_k, 1.0, ux_p, offs, ux_p, offs);
                    // calculate lamda_tilde_k
                    // copy vk to right location
                    // we implemented a version of vector copy that starts with copy of last element, to avoid aliasing error
                    VECCPR(gammamrho_k, lam_p, offs_g_k, lam_p, offs_g_k + rho_k);
                    ROWEX(rho_k, -1.0, RSQrqt_tilde_p + k, nu + nx, 0, lam_p, offs_g_k);
                    // assume aliasing of last two eliments is allowed
                    GEMV_T(nu + nx, rho_k, -1.0, RSQrqt_tilde_p + k, 0, 0, ux_p, offs, 1.0, lam_p, offs_g_k, lam_p, offs_g_k);
                    // nu-rank_k+nx,0
                    // needless copy because feature not implemented yet in trsv_lnn
                    GECP(rho_k, gamma_k, Ggt_tilde_p + k, nu - rho_k + nx + 1, 0, AL_p, 0, 0);
                    // U^-T
                    TRSV_LNN(rho_k, AL_p, 0, 0, lam_p, offs_g_k, lam_p, offs_g_k);
                    // L^-T
                    TRSV_UNU(rho_k, gamma_k, AL_p, 0, 0, lam_p, offs_g_k, lam_p, offs_g_k);
                    (Pl_p + k)->PtV(rho_k, lam_p, offs_g_k);
                    (Pr_p + k)->PtV(rho_k, ux_p, offs);
                }
                // calculate xkp1
                ROWEX(nxp1, 1.0, BAbt_p + k, nu + nx, 0, ux_p, offsp1 + nup1);
                GEMV_T(nu + nx, nxp1, 1.0, BAbt_p + k, 0, 0, ux_p, offs, 1.0, ux_p, offsp1 + nup1, ux_p, offsp1 + nup1);
                // calculate lam_dyn xp1
                ROWEX(nxp1, 1.0, Ppt_p + (k + 1), nxp1, 0, lam_p, offs_dyn_eq_k);
                GEMV_T(nxp1, nxp1, 1.0, Ppt_p + (k + 1), 0, 0, ux_p, offsp1 + nup1, 1.0, lam_p, offs_dyn_eq_k, lam_p, offs_dyn_eq_k);
                GEMV_T(gammamrho_kp1, nxp1, 1.0, Hh_p + (k + 1), 0, 0, lam_p, offs_g_kp1, 1.0, lam_p, offs_dyn_eq_k, lam_p, offs_dyn_eq_k);
                if (ng_ineq > 0)
                {
                    // calculate delta_s
                    ROWEX(ng_ineq, 1.0, Ggt_ineq_p + k, nu + nx, 0, delta_s_p, offs_ineq_k);
                    GEMV_T(nu + nx, ng_ineq, 1.0, Ggt_ineq_p + k, 0, 0, ux_p, offs, 1.0, delta_s_p, offs_ineq_k, delta_s_p, offs_ineq_k);
                    // calculate lamineq
                    for (int i = 0; i < ng_ineq; i++)
                    {
                        double scaling_factor_L = inertia_correction;
                        double scaling_factor_U = inertia_correction;
                        double zLi = VECEL(zL_p, offs_ineq_k + i);
                        double zUi = VECEL(zU_p, offs_ineq_k + i);
                        double si = VECEL(s_p, offs_ineq_k + i);
                        double loweri = VECEL(lower_p, offs_ineq_k + i);
                        double upperi = VECEL(upper_p, offs_ineq_k + i);
                        double grad_barrier_L = 0.0;
                        double grad_barrier_U = 0.0;
                        if (!isinf(loweri))
                        {
                            double dist = si - loweri;
                            double dist_m1 = 1.0 / dist;
                            scaling_factor_L = zLi * dist_m1;
                            grad_barrier_L = mu * dist_m1;
                            VECEL(delta_zL_p, offs_ineq_k + i) = grad_barrier_L - VECEL(zL_p, offs_ineq_k + i) - scaling_factor_L * VECEL(delta_s_p, offs_ineq_k + i);
                        }
                        if (!isinf(upperi))
                        {
                            double dist = upperi - si;
                            double dist_m1 = 1.0 / dist;
                            scaling_factor_U += zUi * dist_m1;
                            grad_barrier_U = -mu * dist_m1;
                            VECEL(delta_zU_p, offs_ineq_k + i) = grad_barrier_U - VECEL(zU_p, offs_ineq_k + i) - scaling_factor_U * VECEL(delta_s_p, offs_ineq_k + i);
                        }
                        VECEL(lam_p, offs_g_ineq_k + i) = grad_barrier_L + grad_barrier_U + (scaling_factor_L + scaling_factor_U) * VECEL(delta_s_p, offs_ineq_k + i);
                    }
                }
            }
            // double el = blasfeo_toc(&timer);
            // cout << "el time " << el << endl;
            return 0;
        }
        FatropMemoryMatBF Ppt;
        FatropMemoryMatBF Hh;
        FatropMemoryMatBF AL;
        FatropMemoryMatBF RSQrqt_tilde;
        FatropMemoryMatBF Ggt_stripe;
        FatropMemoryMatBF Ggt_tilde;
        FatropMemoryMatBF GgLt;
        FatropMemoryMatBF RSQrqt_hat;
        FatropMemoryMatBF Llt;
        FatropMemoryMatBF Llt_shift; // needed because feature not implemented yet
        FatropMemoryMatBF GgIt_tilde;
        FatropMemoryMatBF GgLIt;
        FatropMemoryMatBF HhIt;
        FatropMemoryMatBF PpIt_hat;
        FatropMemoryMatBF LlIt;
        FatropMemoryMatBF Ggt_ineq_temp;
        MemoryPermMat Pl;
        MemoryPermMat Pr;
        MemoryPermMat PlI;
        MemoryPermMat PrI;
        FatropVector<int> gamma;
        FatropVector<int> rho;
    };
};     // namespace
#endif // OCPLSRICCATIINCLUDED