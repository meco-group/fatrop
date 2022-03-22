#ifndef OCPRESTOLSINCLUDED
#define OCPRESTOLSINCLUDED
#include "OCPKKT.hpp"
#include <cmath>
#include "OCPResto.hpp"
namespace fatrop
{
    class OCPRestoLS
    {
    public:
        OCPRestoLS(const OCPDims &dims, const NLPDims &nlpdims) : p_offs_(nlpdims.nineqs),
                                                                  n_offs_(p_offs_ + nlpdims.neqs),
                                                                  primal_s_offs_(nlpdims.neqs),
                                                                  nlpdims_(nlpdims),
                                                                  RSQrqt_tilde(dims.nu + dims.nx + 1, dims.nu + dims.nx, dims.K),
                                                                  BAbt_tmp(dims.nu + dims.nx + 1, rotate(dims.nx, 1), dims.K - 1),
                                                                  Ggt_tmp(dims.nu + dims.nx + 1, dims.ng, dims.K),
                                                                  Ggt_ineq_tmp(dims.nu + dims.nx + 1, dims.ng_ineq, dims.K),
                                                                  L(dims.nu + dims.nx + 1, rotate(dims.nx, 1), dims.K - 1),
                                                                  Ppt(dims.nx + 1, dims.nx, dims.K),
                                                                  ukxkp1(max(dims.nx) + max(dims.nu), 1){};
        int computeSD(
            OCPKKTMemory *OCP,
            const double inertia_correction_w,
            const double inertia_correction_c,
            const double mu,
            const double kappa_d,
            const double rho,
            const double zeta,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &s,
            const FatropVecBF &zL,
            const FatropVecBF &zU,
            const FatropVecBF &delta_zL,
            const FatropVecBF &delta_zU,
            const FatropVecBF &lower,
            const FatropVecBF &upper,
            const FatropVecBF &delta_s,
            const FatropVecBF &grad_xs,
            const FatropVecBF &hess_xs)
        {
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP->name)
#define AUXMACRO(type, name, suffix) type name##suffix = ((type)OCP->aux.name)
#define SOLVERMACRO(type, name, suffix) type name##suffix = ((type)name)
            int K = OCP->K;
            // make variables local for efficiency
            OCPMACRO(MAT *, RSQrqt, _p);
            OCPMACRO(MAT *, BAbt, _p);
            OCPMACRO(MAT *, Ggt, _p);
            OCPMACRO(MAT *, Ggt_ineq, _p);
            OCPMACRO(int *, nu, _p);
            OCPMACRO(int *, nx, _p);
            OCPMACRO(int *, ng, _p);
            OCPMACRO(int *, ng_ineq, _p);
            SOLVERMACRO(MAT *, Ggt_tmp, _p);
            SOLVERMACRO(MAT *, Ggt_ineq_tmp, _p);
            SOLVERMACRO(MAT *, BAbt_tmp, _p);
            SOLVERMACRO(MAT *, RSQrqt_tilde, _p);
            SOLVERMACRO(MAT *, L, _p);
            SOLVERMACRO(MAT *, Ppt, _p);
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
            SOLVERMACRO(VEC *, grad_xs, _p);
            SOLVERMACRO(VEC *, hess_xs, _p);
            SOLVERMACRO(VEC *, ukxkp1, _p);
            int *offs_ineq_p = (int *)OCP->aux.ineq_offs.data();
            int *offs_ux_p = (int *)OCP->aux.ux_offs.data();
            int *offs_g_p = (int *)OCP->aux.g_offs.data();
            int *offs_dyn_eq_p = (int *)OCP->aux.dyn_eq_offs.data();
            int *offs_g_ineq_p = (int *)OCP->aux.g_ineq_offs.data();
            const int p_offs = p_offs_;
            const int n_offs = n_offs_;
            const int n_nlp_vars = nlpdims_.nvars;
            const int n_nlp_ineqs = nlpdims_.nineqs;
            const int n_nlp_eqs = nlpdims_.neqs;

///////////////////////////////////////////////////////////////////////////////////////
//// note we calculate dlam instead of lam, like in the ipopt implementation paper ////
///////////////////////////////////////////////////////////////////////////////////////
#define SIGMA_MACRO double sigma_pnm1 = 1.0 / (pi / zpi + ni / zni)

            // last stage
            {
                const int nx = nx_p[K - 1];
                const int nu = nu_p[K - 1]; // this should be zero but is included here in case of misuse
                const int ng = ng_p[K - 1];
                const int ng_ineq = ng_ineq_p[K - 1];
                const int offs_ineq_k = offs_ineq_p[K - 1];
                const int offs_g_ineq_k = offs_g_ineq_p[K - 1];
                const int offs_ux_k = offs_ux_p[K - 1];
                const int offs_g_k = offs_g_p[K - 1];
                GECP(nx + 1, nx, RSQrqt_p + K - 1, nu, nu, RSQrqt_tilde_p + K - 1, 0, 0);
                // RSQrqt
                DIAAD(nx, 1.0, hess_xs_p, offs_ux_k + nu, RSQrqt_tilde_p + K - 1, 0, 0);
                for (int i = 0; i < nx; i++)
                {
                    const double grad_i = VECEL(grad_xs_p, offs_ux_k + nu + i);
                    MATEL(RSQrqt_tilde_p + K - 1, nu + nx, i + nu) = grad_i;
                }
                // Ggt
                GECP(nx + 1, ng, Ggt_p + K - 1, nu, 0, Ggt_tmp_p + K - 1, 0, 0);
                for (int i = 0; i < ng; i++)
                {
                    const double pi = VECEL(s_p, p_offs + offs_g_k + i);
                    const double ni = VECEL(s_p, n_offs + offs_g_k + i);
                    const double zpi = VECEL(zL_p, p_offs + offs_g_k + i);
                    const double zni = VECEL(zL_p, n_offs + offs_g_k + i);
                    // double tmp = -pi + ni + (rho * (mu - pi) + pi) / zpi + (rho * (mu - ni) + ni) / zni;
                    double tmp = -pi + ni + (rho * (mu - pi) ) / zpi + (rho * (mu - ni) ) / zni;
                    // MATEL(Ggt_p + K - 1, nu + nx, i) += tmp;
                    MATEL(Ggt_tmp_p + K - 1, nu + nx, i) += tmp;
                    SIGMA_MACRO;
                    // prepare scaled version of Ggt
                    COLSC(nx + 1, sigma_pnm1, Ggt_tmp_p + K - 1, 0, i);
                }
                // add the penalty
                SYRK_LN_MN(nx + 1, nx, ng, 1.0, Ggt_tmp_p + K - 1, 0, 0, Ggt_p + K - 1, nu, 0, 1.0, RSQrqt_tilde_p + K - 1, 0, 0, RSQrqt_tilde_p + K - 1, 0, 0);
                // Ggt_ineq
                if (ng_ineq > 0)
                {
                    GECP(nx + 1, ng_ineq, Ggt_ineq_p + K - 1, nu, 0, Ggt_ineq_tmp_p + K - 1, 0, 0);
                    for (int i = 0; i < ng_ineq; i++)
                    {
                        // calculation of grad_barrier and scaling factor
                        double scaling_factor = inertia_correction_w;
                        const double zLi = VECEL(zL_p, offs_ineq_k + i);
                        const double zUi = VECEL(zU_p, offs_ineq_k + i);
                        const double si = VECEL(s_p, offs_ineq_k + i);
                        const double loweri = VECEL(lower_p, offs_ineq_k + i);
                        const double upperi = VECEL(upper_p, offs_ineq_k + i);
                        bool lower_bounded = !isinf(loweri);
                        bool upper_bounded = !isinf(upperi);
                        double grad_barrier = 0.0;
                        if (lower_bounded)
                        {
                            double dist = si - loweri;
                            double dist_m1 = 1.0 / dist;
                            scaling_factor += zLi * dist_m1;
                            grad_barrier -= mu * dist_m1;
                        }
                        if (upper_bounded)
                        {
                            double dist = upperi - si;
                            double dist_m1 = 1.0 / dist;
                            scaling_factor += zUi * dist_m1;
                            grad_barrier += mu * dist_m1;
                        }
                        if (!(lower_bounded && upper_bounded))
                        {
                            grad_barrier += lower_bounded ? kappa_d * mu : -kappa_d * mu;
                        }
                        // calculation of gi
                        const double pi = VECEL(s_p, p_offs + offs_g_ineq_k + i);
                        const double ni = VECEL(s_p, n_offs + offs_g_ineq_k + i);
                        const double zpi = VECEL(zL_p, p_offs + offs_g_ineq_k + i);
                        const double zni = VECEL(zL_p, n_offs + offs_g_ineq_k + i);
                        // double tmp = -pi + ni + (rho * (mu - pi) + pi) / zpi + (rho * (mu - ni) + ni) / zni;
                    double tmp = -pi + ni + (rho * (mu - pi) ) / zpi + (rho * (mu - ni) ) / zni;
                        double gi = MATEL(Ggt_ineq_p + K - 1, nu + nx, i) + tmp;
                        double hess_i = VECEL(hess_xs_p, n_nlp_vars + offs_ineq_k + i);
                        double sigma_pnm1 = 1.0 / (pi / zpi + ni / zni);
                        double sigma_stripe_i = scaling_factor + sigma_pnm1 + hess_i;
                        tmp = sigma_pnm1 * sigma_pnm1 / sigma_stripe_i;
                        double sigma_a = sigma_pnm1 - tmp;
                        double sigma_b = sigma_pnm1 + tmp;
                        COLSC(nx, sigma_a, Ggt_ineq_tmp_p + K - 1, 0, i);
                        MATEL(Ggt_ineq_tmp_p + K - 1, nx, i) = (grad_barrier - sigma_pnm1 * gi) * sigma_b + gi * sigma_pnm1;
                    }
                    // add the penalty
                    // SYRK_LN_MN(nx + 1, nx, ng_ineq, 1.0, Ggt_ineq_tmp_p + K - 1, 0, 0, Ggt_ineq_p + K - 1, nu, 0, 1.0, RSQrqt_tilde_p + K - 1, 0, 0, RSQrqt_tilde_p + K - 1, 0, 0);
                }
            }
            // other stages
            for (int k = K - 2; k >= 0; --k)
            {
                const int nu = nu_p[k];
                const int nup1 = nu_p[k + 1];
                const int nx = nx_p[k];
                const int nxp1 = nx_p[k + 1];
                const int ng = ng_p[k];
                const int ng_ineq = ng_ineq_p[k];
                const int offs_ineq_k = offs_ineq_p[k];
                const int offs_g_ineq_k = offs_g_ineq_p[k];
                const int offs_ux_k = offs_ux_p[k];
                const int offs_g_k = offs_g_p[k];
                const int offs_dyn_eq_k = offs_dyn_eq_p[k];
                // RSQrqt
                GECP(nu + nx + 1, nu + nx, RSQrqt_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
                DIAAD(nu + nx, 1.0, hess_xs_p, offs_ux_k, RSQrqt_tilde_p + k, 0, 0);
                for (int i = 0; i < nu + nx; i++)
                {
                    const double grad_i = VECEL(grad_xs_p, offs_ux_k + i);
                    MATEL(RSQrqt_tilde_p + k, nu + nx, i) = grad_i;
                }
                // Ggt
                GECP(nu + nx + 1, ng, Ggt_p + k, 0, 0, Ggt_tmp_p + k, 0, 0);
                for (int i = 0; i < ng; i++)
                {
                    const double pi = VECEL(s_p, p_offs + offs_g_k + i);
                    const double ni = VECEL(s_p, n_offs + offs_g_k + i);
                    const double zpi = VECEL(zL_p, p_offs + offs_g_k + i);
                    const double zni = VECEL(zL_p, n_offs + offs_g_k + i);
                    // double tmp = -pi + ni + (rho * (mu - pi) + pi) / zpi + (rho * (mu - ni) + ni) / zni;
                    double tmp = -pi + ni + (rho * (mu - pi) ) / zpi + (rho * (mu - ni) ) / zni;
                    // MATEL(Ggt_p + k, nu + nx, i) += tmp;
                    MATEL(Ggt_tmp_p + k, nu + nx, i) += tmp;
                    // prepare scaled version of Ggt
                    SIGMA_MACRO;
                    COLSC(nu + nx + 1, sigma_pnm1, Ggt_tmp_p + k, 0, i);
                }
                SYRK_LN_MN(nu + nx + 1, nu + nx, ng, 1.0, Ggt_tmp_p + k, 0, 0, Ggt_p + k, 0, 0, 1.0, RSQrqt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
                // Ggt_ineq
                if (ng_ineq > 0)
                {
                    GECP(nu + nx + 1, ng_ineq, Ggt_ineq_p + k, 0, 0, Ggt_ineq_tmp_p + k, 0, 0);
                    for (int i = 0; i < ng_ineq; i++)
                    {
                        // calculation of grad_barrier and scaling factor
                        double scaling_factor = inertia_correction_w;
                        const double zLi = VECEL(zL_p, offs_ineq_k + i);
                        const double zUi = VECEL(zU_p, offs_ineq_k + i);
                        const double si = VECEL(s_p, offs_ineq_k + i);
                        const double loweri = VECEL(lower_p, offs_ineq_k + i);
                        const double upperi = VECEL(upper_p, offs_ineq_k + i);
                        double grad_barrier = 0.0;
                        bool lower_bounded = !isinf(loweri);
                        bool upper_bounded = !isinf(upperi);
                        if (lower_bounded)
                        {
                            double dist = si - loweri;
                            double dist_m1 = 1.0 / dist;
                            scaling_factor += zLi * dist_m1;
                            grad_barrier -= mu * dist_m1;
                        }
                        if (upper_bounded)
                        {
                            double dist = upperi - si;
                            double dist_m1 = 1.0 / dist;
                            scaling_factor += zUi * dist_m1;
                            grad_barrier += mu * dist_m1;
                        }
                        if (!(lower_bounded && upper_bounded))
                        {
                            grad_barrier += lower_bounded ? kappa_d * mu : -kappa_d * mu;
                        }
                        // calculation of gi
                        const double pi = VECEL(s_p, p_offs + offs_g_ineq_k + i);
                        const double ni = VECEL(s_p, n_offs + offs_g_ineq_k + i);
                        const double zpi = VECEL(zL_p, p_offs + offs_g_ineq_k + i);
                        const double zni = VECEL(zL_p, n_offs + offs_g_ineq_k + i);
                        // double tmp = -pi + ni + (rho * (mu - pi) + pi) / zpi + (rho * (mu - ni) + ni) / zni;
                        double tmp = -pi + ni + (rho * (mu - pi) ) / zpi + (rho * (mu - ni)) / zni;
                        double gi = MATEL(Ggt_ineq_p + k, nu + nx, i) + tmp;
                        double hess_i = VECEL(hess_xs_p, n_nlp_vars + offs_ineq_k + i);
                        double sigma_pnm1 = 1.0 / (pi / zpi + ni / zni);
                        double sigma_stripe_i = scaling_factor + sigma_pnm1 + hess_i;
                        tmp = sigma_pnm1 * sigma_pnm1 / sigma_stripe_i;
                        double sigma_a = sigma_pnm1 - tmp;
                        double sigma_b = sigma_pnm1 + tmp;
                        COLSC(nu + nx, sigma_a, Ggt_ineq_tmp_p + k, 0, i);
                        MATEL(Ggt_ineq_tmp_p + k, nu + nx, i) = (grad_barrier - sigma_pnm1 * gi) * sigma_b + gi * sigma_pnm1;
                    }
                    // add the penalty
                    // SYRK_LN_MN(nu + nx + 1, nu+nx, ng_ineq, 1.0, Ggt_ineq_tmp_p + k, 0, 0, Ggt_ineq_p + k, 0, 0, 1.0, RSQrqt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
                }
                GECP(nu + nx + 1, nxp1, BAbt_p + k, 0, 0, BAbt_tmp_p + k, 0, 0);
                // dynamics eqs
                for (int i = 0; i < nxp1; i++)
                {
                    const double pi = VECEL(s_p, p_offs + offs_dyn_eq_k + i);
                    const double ni = VECEL(s_p, n_offs + offs_dyn_eq_k + i);
                    const double zpi = VECEL(zL_p, p_offs + offs_dyn_eq_k + i);
                    const double zni = VECEL(zL_p, n_offs + offs_dyn_eq_k + i);
                    // double tmp = -pi + ni + (rho * (mu - pi) + pi) / zpi + (rho * (mu - ni) + ni) / zni;
                    double tmp = -pi + ni + (rho * (mu - pi) ) / zpi + (rho * (mu - ni) ) / zni;
                    // MATEL(BAbt_p + k, nu + nx, i) += tmp;
                    MATEL(BAbt_tmp_p + k, nu + nx, i) += tmp;
                    SIGMA_MACRO;
                    // [-I B A] --> -I
                    MATEL(RSQrqt_tilde_p + k + 1, nup1 + i, nup1 + i) += sigma_pnm1;
                    MATEL(RSQrqt_tilde_p + k + 1, nup1 + nxp1, nup1 + i) -= sigma_pnm1 * MATEL(BAbt_p + k, nu + nx, i);
                    // prepare scaled version of BAbt
                    COLSC(nu + nx + 1, sigma_pnm1, BAbt_tmp_p + k, 0, i);
                }
                // add the penalty
                SYRK_LN_MN(nu + nx + 1, nu + nx, nxp1, 1.0, BAbt_tmp_p + k, 0, 0, BAbt_p + k, 0, 0, 1.0, RSQrqt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
            }
            {
                // save Pk
                const int nu = nu_p[K - 1];
                const int nx = nx_p[K - 1];
                GECP(nx + 1, nx, RSQrqt_tilde_p + K - 1, nu, nu, Ppt_p + K - 1, 0, 0);
                POTRF_L_MN(nx + 1, nx, Ppt_p + K - 1, 0, 0, Ppt_p + K - 1, 0, 0);
            }
            // calculation of delta x -> backward recursion
            for (int k = K - 2; k >= 0; --k)
            {
                const int nu = nu_p[k];
                const int nx = nx_p[k];
                const int nxp1 = nx_p[k + 1];

                /* Lmatrix structure
                -BTS
                -ATS
                -beTS
                ------
                Pk+1 = LLT
                ----->
                RSQrat_hat_k - A^T L-T L-1 A (SYRK)
                */
                GECPSC(nu + nx + 1, nxp1, -1.0, BAbt_p + k, 0, 0, L_p + k, 0, 0);
                // assume aliasing is allowed for last two elements
                TRSM_RLTN(nu + nx + 1, nxp1, 1.0, Ppt_p + k + 1, 0, 0, L_p + k, 0, 0, L_p + k, 0, 0);
                // update RSQrqt_hat
                SYRK_LN_MN(nu + nx + 1, nu + nx, nxp1, -1.0, L_p + k, 0, 0, L_p + k, 0, 0, 1.0, RSQrqt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
                GECP(nx + 1, nx, RSQrqt_tilde_p + k, nu, nu, Ppt_p + k, 0, 0);
            }

            // u0x0
            {
                const int nu = nu_p[0];
                const int nx = nx_p[0];
                ROWEX(nx + nu, -1.0, RSQrqt_tilde_p + 0, nu + nx, 0, ux_p, 0);
            }
            for (int k = 0; k < K - 1; k++)
            {
                // uk xk+1
                const int nu = nu_p[k];
                const int nup1 = nu_p[k + 1];
                const int nx = nx_p[k];
                const int nxp1 = nx_p[k + 1];
                const int offs_ux_k = offs_ux_p[k];
                const int offs_ux_kp1 = offs_ux_p[k + 1];
                TRSV_LTN(nx + nu, RSQrqt_tilde_p + k, 0, 0, ux_p, offs_ux_k, ux_p, offs_ux_k);
                GEMV_T(nu + nx, nxp1, -1.0, L_p + k, 0, 0, ux_p, offs_ux_k + nu, 1.0, ux_p, offs_ux_kp1 + nup1, ux_p, offs_ux_kp1 + nup1);
                ROWEX(nxp1 + nup1, -1.0, RSQrqt_tilde_p + k + 1, nup1 + nxp1, 0, ux_p, offs_ux_kp1);
            }
            {
                const int nu = nu_p[K - 1];
                const int nx = nx_p[K - 1];
                const int offs_ux_k = offs_ux_p[K - 1];
                TRSV_LTN(nx ,RSQrqt_tilde_p + K - 1, nu, nu, ux_p, offs_ux_k + nu, ux_p, offs_ux_k + nu);
            }

            // calculate lam_ineq + deltas
            for (int k = 0; k < K - 1; k++)
            {
                const int nu = nu_p[k];
                const int nx = nx_p[k];
                const int ng_ineq = ng_ineq_p[k];
                const int offs_ineq_k = offs_ineq_p[k];
                const int offs_g_ineq_k = offs_g_ineq_p[k];
                const int offs_ux_k = offs_ux_p[k];
                if (ng_ineq > 0)
                {
                    ROWEX(ng_ineq, 1.0, Ggt_ineq_p, nu + nx, 0, delta_s_p, offs_ineq_k);
                    GEMV_T(nu + nx + 1, ng_ineq, 1.0, Ggt_ineq_p + k, 0, 0, ux_p, offs_ux_k, 1.0, delta_s_p, offs_ineq_k, delta_s_p, offs_ineq_k);
                    for (int i = 0; i < ng_ineq; i++)
                    {
                        // calculation of grad_barrier and scaling factor
                        double scaling_factor = inertia_correction_w;
                        const double zLi = VECEL(zL_p, offs_ineq_k + i);
                        const double zUi = VECEL(zU_p, offs_ineq_k + i);
                        const double si = VECEL(s_p, offs_ineq_k + i);
                        const double loweri = VECEL(lower_p, offs_ineq_k + i);
                        const double upperi = VECEL(upper_p, offs_ineq_k + i);
                        double grad_barrier = 0.0;
                        if (!isinf(loweri))
                        {
                            double dist = si - loweri;
                            double dist_m1 = 1.0 / dist;
                            scaling_factor += zLi * dist_m1;
                            grad_barrier -= mu * dist_m1;
                        }
                        if (!isinf(upperi))
                        {
                            double dist = upperi - si;
                            double dist_m1 = 1.0 / dist;
                            scaling_factor += zUi * dist_m1;
                            grad_barrier += mu * dist_m1;
                        }
                        // calculation of gi
                        const double pi = VECEL(s_p, p_offs + offs_g_ineq_k + i);
                        const double ni = VECEL(s_p, n_offs + offs_g_ineq_k + i);
                        const double zpi = VECEL(zL_p, p_offs + offs_g_ineq_k + i);
                        const double zni = VECEL(zL_p, n_offs + offs_g_ineq_k + i);
                        double hess_i = VECEL(hess_xs_p, n_nlp_vars + offs_ineq_k + i);
                        double sigma_pnm1 = 1.0 / (pi / zpi + ni / zni);
                        double sigma_stripe_i = scaling_factor + sigma_pnm1 + hess_i;
                        double tmp2 = sigma_pnm1 * VECEL(delta_s_p, offs_ineq_k + i);
                        double dsi = -1.0 / sigma_stripe_i * (-tmp2 + grad_barrier);
                        VECEL(delta_s_p, offs_ineq_k + i) = dsi;
                        VECEL(lam_p, offs_g_ineq_k + i) = tmp2 - sigma_pnm1 * dsi;
                    }
                }
            }
            // lam_ineq + delta_s last stage
            {
                const int nu = nu_p[K - 1];
                const int nx = nx_p[K - 1];
                const int ng_ineq = ng_ineq_p[K - 1];
                const int offs_ineq_k = offs_ineq_p[K - 1];
                const int offs_g_ineq_k = offs_g_ineq_p[K - 1];
                const int offs_ux_k = offs_ux_p[K - 1];
                if (ng_ineq > 0)
                {
                    ROWEX(ng_ineq, 1.0, Ggt_ineq_p, nu + nx, 0, delta_s_p, offs_ineq_k);
                    GEMV_T(nx + 1, ng_ineq, 1.0, Ggt_ineq_p + K - 1, nu, 0, ux_p, offs_ux_k + nu, 1.0, delta_s_p, offs_ineq_k, delta_s_p, offs_ineq_k);
                    for (int i = 0; i < ng_ineq; i++)
                    {
                        // calculation of grad_barrier and scaling factor
                        double scaling_factor = inertia_correction_w;
                        const double zLi = VECEL(zL_p, offs_ineq_k + i);
                        const double zUi = VECEL(zU_p, offs_ineq_k + i);
                        const double si = VECEL(s_p, offs_ineq_k + i);
                        const double loweri = VECEL(lower_p, offs_ineq_k + i);
                        const double upperi = VECEL(upper_p, offs_ineq_k + i);
                        bool lower_bounded = !isinf(loweri);
                        bool upper_bounded = !isinf(upperi);
                        double grad_barrier = 0.0;
                        if (lower_bounded)
                        {
                            double dist = si - loweri;
                            double dist_m1 = 1.0 / dist;
                            scaling_factor += zLi * dist_m1;
                            grad_barrier -= mu * dist_m1;
                        }
                        if (upper_bounded)
                        {
                            double dist = upperi - si;
                            double dist_m1 = 1.0 / dist;
                            scaling_factor += zUi * dist_m1;
                            grad_barrier += mu * dist_m1;
                        }
                        if (!(lower_bounded && upper_bounded))
                        {
                            grad_barrier += lower_bounded ? kappa_d * mu : -kappa_d * mu;
                        }
                        // calculation of gi
                        const double pi = VECEL(s_p, p_offs + offs_g_ineq_k + i);
                        const double ni = VECEL(s_p, n_offs + offs_g_ineq_k + i);
                        const double zpi = VECEL(zL_p, p_offs + offs_g_ineq_k + i);
                        const double zni = VECEL(zL_p, n_offs + offs_g_ineq_k + i);
                        double hess_i = VECEL(hess_xs_p, n_nlp_vars + offs_ineq_k + i);
                        double sigma_pnm1 = 1.0 / (pi / zpi + ni / zni);
                        double sigma_stripe_i = scaling_factor + sigma_pnm1 + hess_i;
                        double tmp2 = sigma_pnm1 * VECEL(delta_s_p, offs_ineq_k + i);
                        double dsi = -1.0 / sigma_stripe_i * (-tmp2 + grad_barrier);
                        VECEL(delta_s_p, offs_ineq_k + i) = dsi;
                        VECEL(lam_p, offs_g_ineq_k + i) = tmp2 - sigma_pnm1 * dsi;
                    }
                }
            }
            // calculate lams of stagewise constraints
            for (int k = 0; k < K; k++)
            {
                const int nu = nu_p[k];
                const int nx = nx_p[k];
                const int ng = ng_p[k];
                const int offs_ux_k = offs_ux_p[k];
                const int offs_g_k = offs_g_p[k];
                ROWEX(ng, 1.0, Ggt_tmp_p + k, 0, 0, lam_p, offs_g_k);
                GEMV_T(nu + nx, ng, 1.0, Ggt_tmp_p + k, 0, 0, ux_p, offs_ux_k, 1.0, lam_p, offs_g_k, lam_p, offs_g_k);
            }
            // calculate lams of dynamics constraints
            for (int k = 0; k < K - 1; k++)
            {
                const int nu = nu_p[k];
                const int nx = nx_p[k];
                const int nxp1 = nx_p[k + 1];
                const int offs_ux_k = offs_ux_p[k];
                const int offs_ux_kp1 = offs_ux_p[k + 1];
                const int offs_dyn_eq_k = offs_dyn_eq_p[k];
                ROWEX(nxp1, 1.0, BAbt_tmp_p + k, 0, 0, lam_p, offs_dyn_eq_k);
                AXPY(nxp1, -1.0, ux_p, offs_ux_kp1, lam_p, offs_dyn_eq_k, lam_p, offs_dyn_eq_k);
                GEMV_T(nu + nx, nxp1, 1.0, BAbt_tmp_p + k, 0, 0, ux_p, offs_ux_k, 1.0, lam_p, offs_dyn_eq_k, lam_p, offs_dyn_eq_k);
            }
            // calculate delta p, delta n, delta zp, delta zn
            for (int i = 0; i < n_nlp_eqs; i++)
            {
                const double lami = VECEL(lam_p, i);
                const double pi = VECEL(s_p, p_offs + i);
                const double ni = VECEL(s_p, n_offs + i);
                const double zpi = VECEL(zL_p, p_offs + i);
                const double zni = VECEL(zL_p, n_offs + i);
                const double dp = (mu + pi * (lami - rho)) / zpi;
                const double dn = (mu + ni * (lami - rho)) / zni;
                const double dzp = mu / pi - zpi - zpi / pi * dp;
                const double dzn = mu / ni - zni - zni / ni * dn;
                VECEL(delta_s_p, p_offs + i) = dp;
                VECEL(delta_s_p, n_offs + i) = dn;
                VECEL(delta_zL_p, p_offs + i) = dzp;
                VECEL(delta_zL_p, n_offs + i) = dzn;
            }
            // calculate delta zs
            for (int i = 0; i < n_nlp_ineqs; i++)
            {
                double loweri = VECEL(lower_p, +i);
                double upperi = VECEL(upper_p, +i);
                bool lower_bounded = !isinf(loweri);
                bool upper_bounded = !isinf(upperi);
                double ds = VECEL(delta_s_p, i);
                double zLi = VECEL(zL_p, i);
                double zUi = VECEL(zU_p, i);
                double si = VECEL(s_p, i);
                if (lower_bounded)
                {
                    double dist = si - loweri;
                    double dist_m1 = 1.0 / dist;
                    double scaling_factor_L = zLi * dist_m1;
                    double grad_barrier_L = -mu * dist_m1;
                    double z = VECEL(zL_p, i);
                    double dz = -grad_barrier_L - z - scaling_factor_L * ds;
                    VECEL(delta_zL_p, i) = dz;
                }
                if (upper_bounded)
                {
                    double dist = upperi - si;
                    double dist_m1 = 1.0 / dist;
                    double scaling_factor_U = zUi * dist_m1;
                    double grad_barrier_U = mu * dist_m1;
                    double z = VECEL(zU_p, i);
                    double dz = grad_barrier_U - z + scaling_factor_U * ds;
                    VECEL(delta_zU_p, i) = dz;
                }
            }
            // copy slack variables to primal variables
            VECCP(n_nlp_ineqs + 2 * n_nlp_eqs, delta_s_p, 0, ux_p, n_nlp_vars);
            cout << "ux \n";
            ux.print();
            // cout << "zL \n";
            // zL.print();
            return 0;
        }
        int p_offs_;
        int n_offs_;
        int primal_s_offs_;
        NLPDims nlpdims_;
        FatropMemoryMatBF RSQrqt_tilde;
        FatropMemoryMatBF BAbt_tmp;
        FatropMemoryMatBF Ggt_tmp;
        FatropMemoryMatBF Ggt_ineq_tmp;
        FatropMemoryMatBF L;
        FatropMemoryMatBF Ppt;
        FatropMemoryVecBF ukxkp1;
    };
} // namespace fatrop
#endif // OCPRESTOLSINCLUDED