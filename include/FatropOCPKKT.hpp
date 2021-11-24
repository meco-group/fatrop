/**
 * @file FatropOCPKKT.hpp
 * @author your name (you@domain.com)
 * @brief this file contains classes to represent a KKT matrix by blasfeo submatrices 
 * @version 0.1
 * @date 2021-11-10
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef FATROPOCPKKTINDLUCED
#define FATROPOCPKKTINDLUCED
#include "FatropLinearAlgebraBlasfeo.hpp"
#include "FatropOCP.hpp"
#include "FatropAux.hpp"
using namespace std;
namespace fatrop
{
    /** \brief this class contains all information to represent the KKT system of an equality constrained OCP*/
    class OCP_KKT
    {
    public:
        OCP_KKT(const OCP_dims &dims, fatrop_memory_allocator &fma) : aux(dims, fma), K(dims.K), nu(dims.K, vector<int>(dims.nu), fma), nx(dims.K, vector<int>(dims.nx), fma), ng(dims.K, vector<int>(dims.ng), fma), RSQrqt(vector<int>(dims.nu + dims.nx + 1), vector<int>(dims.nu + dims.nx), dims.K, fma), BAbt(vector<int>(dims.nu + dims.nx + 1), vector<int>(dims.nx), dims.K, fma), Ggt(vector<int>(dims.nu + dims.nx + 1), vector<int>(dims.ng), dims.K, fma){};
        int K;
        fatrop_memory_el<int> nu;
        fatrop_memory_el<int> nx;
        fatrop_memory_el<int> ng;
        /// small-scale Hessian
        fatrop_memory_matrix_bf RSQrqt;
        /// small-scale Jacobian dynamics
        fatrop_memory_matrix_bf BAbt;
        /// small-scale Jacobian stagewise eq constraints
        fatrop_memory_matrix_bf Ggt;
        class OCP_aux
        {
        public:
            OCP_aux(const OCP_dims &dims, fatrop_memory_allocator &fma) : ux_offs(dims.K, offsets(dims.nx + dims.nu), fma), g_offs(dims.K, offsets(dims.ng), fma), max_nu(max(dims.nu)), max_nx(max(dims.nx)), max_ng(max(dims.ng)){};
            /// offset arrays are used for efficiency
            fatrop_memory_el<int> ux_offs;
            /// offset arrays are used for efficiency
            fatrop_memory_el<int> g_offs;
            int max_nu;
            int max_nx;
            int max_ng;
        };
        OCP_aux aux;
    };
    class OCP_KKT_solver
    {
    public:
        OCP_KKT_solver(const OCP_dims &dims, fatrop_memory_allocator &fma) : Ppt(dims.nx + 1, vector<int>(dims.nx), dims.K, fma),
                                                                             Hht(vector<int>(1,max(dims.nu + dims.nx + 1)), vector<int>(?????), 1, fma),
                                                                             gamma(dims.K, vector<int>(dims.K, 0), fma),
                                                                             rho(dims.K, vector<int>(dims.K, 0), fma){};
        // solve a KKT system
        void fact_solve(OCP_KKT *OCP, VEC *ux, VEC *lambda)
        {
            // define compiler macros for notational convenience
#define OCPMACRO(type, name, subfix) type name##subfix = ((type)OCP->name)
#define AUXMACRO(type, name, subfix) type name##subfix = ((type)OCP->aux.name)
#define SOLVERMACRO(type, name, subfix) type name##subfix = ((type)name)
            int K = OCP->K;
            // make variables local for efficiency
            OCPMACRO(MAT *, RSQrqt, _p);
            OCPMACRO(MAT *, BAbt, _p);
            OCPMACRO(MAT *, Ggt, _p);
            SOLVERMACRO(MAT *, Ppt, _p);
            SOLVERMACRO(MAT *, Hht, _p);
            AUXMACRO(int, max_nu, );
            AUXMACRO(int, max_nx, );
            AUXMACRO(int, max_ng, );
            AUXMACRO(int *, ux_offs, _p);
            AUXMACRO(int *, g_offs, _p);
            OCPMACRO(int *, nu, _p);
            OCPMACRO(int *, nx, _p);
            OCPMACRO(int *, ng, _p);
            // recursion
            {
                // last stage
                int nx = nx_p[K - 1];
                int nu = nu_p[K - 1]; // this should be zero but is included here in case of misuse
                int ng = ng_p[K - 1]; // this should be zero but is included here in case of misuse
                // Pp_Km1 <- Qq_Km1
                GECP(nx + 1, nx, RSQrqt_p + (K - 1), nu, nu, Ppt_p + K - 1, 0, 0);
                // Hh_Km1 <- Gg_Km1
                GECP(nx + 1, nx, RSQrqt_p + (K - 1), nu, nu, Hht_p, 0, 0);
            }
            for (int k = K - 2; k >= 0; --k)
            {
                int nu = nu_p[k];
                int nx = nx_p[k];
                int ng = ng_p[k];
                goto SUBSDYN;
            SUBSDYN:
                // AL <- [BAb]^T_k P_kp1
                // AL[-1,:] <- AL[-1,:] + p_kp1^T
                // RSQrqt_stripe <- AL[BA] + RSQrqt
                // if ng[k]>0
                // Ggt_stripe  <- Ggt_k
                // if Hkp1 nonempty
                // Ggt_stripe <- [Ggt_k [BAb_k^T]H_kp1]
                // Ggt_stripe[-1,:] <- Ggt_stripe[-1,:] + h_kp1^T
                // gamma_k <- number of eqs represented by Ggt_stripe
            TRANSFORM_AND_SUBSEQ:
                // symmetric transformation, done a little different than in paper, in order to fuse LA operations
                // LU_FACT_TRANSPOSE(Ggtstripe[:gamma_k, nu+nx+1], nu max)
                // Ggt_tilde_k <- Ggt_stripe[rho_k:nu+nx+1, :rho] L-T (note that this is slightly different from the implementation)
                // permutations
                // Hh_k <- Ggt_stripe[nu:nu+nx+1, rho:] (transfer to next stage)
                // GL <- Ggt_tilde_k @ RSQrqt[:, :rho]^T + RSQrqt[rho:nu+nx+1, rho:]^T (note the transpose of the last term!!)
                // RSQrqt_hat = Gt_tilde_k[:-1,:] @ GL[:, :rho]^T + GL[:rho:]^T (note the transpose of the last term!!)
            SCHUR:
                // DLlt_k = [chol(R_hatk)  Llk@chol(R_hatk)^-T]
                // Pp_k = Qq_hatk - L_k^T @ Ll_k
            FIRST_STAGE:
                // if gamma_0 - rho_0 > 0
                // LU_FACT_TRANSPOSE(Hh0)
                // permutations
                // Ggt_tilde_I <- Ggt_stripe[rho_I:nx+1,:rho_I] L^-T
                // h_tilde_I <- - U_I ^-1 Ggt_tilde_I[nx+1, :]
                // ?? if nx - rho_I > 0 ??
                // GL_I <- Ggt_tildeI  @ Ppt_I^T + Ppt[rho:nx+1, rho:]^T
                // Pphat_I <- Ggt_tilde_I[:-1,:] @  GL_I[:,:rho]^T + GL[:rho,:]^T
                // DlI = [chol(Phat_I) lI@chol(phat_I)^-T]
            FORWARD_SUBSTITUTION:
                cout << "test" << endl;
            }
        }
        fatrop_memory_matrix_bf Ppt;
        fatrop_memory_matrix_bf Hht;
        fatrop_memory_el<int> gamma;
        fatrop_memory_el<int> rho;
    };
} // namespace fatrop
#endif // FATROPOCPKKTINDLUCED