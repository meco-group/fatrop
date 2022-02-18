#ifndef DUINFEVALINCLUDED
#define DUINFEVALINCLUDED
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "OCPKKT.hpp"
#include "aux/Common.hpp"
namespace fatrop
{
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP->name)
    class DuInfEvaluator
    {
    public:
        int DuInfEval(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &lam,
            const FatropVecBF &grad,
            FatropVecBF &du_inf)
        {
            int K = OCP->K;
            VEC *lam_p = (VEC *)lam;
            VEC *grad_p = (VEC *)grad;
            VEC *du_inf_p = (VEC *)du_inf;
            OCPMACRO(MAT *, BAbt, _p);
            OCPMACRO(MAT *, Ggt, _p);
            OCPMACRO(MAT *, Ggt_ineq, _p);
            OCPMACRO(int *, nu, _p);
            OCPMACRO(int *, nx, _p);
            OCPMACRO(int *, ng, _p);
            OCPMACRO(int *, ng_ineq, _p);
            DBGASSERT(grad_p->m == du_inf_p->m);
            // copy grad_f to du_inf
            VECCP(grad_p->m, grad_p, 0, du_inf_p, 0);
            int *offs_ux = (int *)OCP->aux.ux_offs.data();
            int *offs_g = (int *)OCP->aux.g_offs.data();
            int *offs_dyn_eq = (int *)OCP->aux.dyn_eq_offs.data();
            int *offs_g_ineq = (int *)OCP->aux.g_ineq_offs.data();
            // contribution of dynamics constraints
            for (int k = 0; k < K - 1; k++)
            {
                const int nu = nu_p[k];
                const int nup1 = nu_p[k + 1];
                const int nx = nx_p[k];
                const int nxp1 = nx_p[k + 1];
                const int offsp1 = offs_ux[k + 1];
                const int offs = offs_ux[k];
                const int offs_dyn_eq_k = offs_dyn_eq[k];
                GEMV_N(nu + nx, nxp1, 1.0, BAbt_p + k, 0, 0, lam_p, offs_dyn_eq_k, 1.0, du_inf_p, offs, du_inf_p, offs);
                AXPY(nxp1, -1.0, lam_p, offs_dyn_eq_k, du_inf_p, offsp1 + nup1, du_inf_p, offsp1 + nup1);
            }
            // contribution of equality constraints
            for (int k = 0; k < K; k++)
            {
                const int nu = nu_p[k];
                const int nx = nx_p[k];
                const int ng = ng_p[k];
                const int offs = offs_ux[k];
                const int offs_g_k = offs_g[k];
                GEMV_N(nu + nx, ng, 1.0, Ggt_p + k, 0, 0, lam_p, offs_g_k, 1.0, du_inf_p, offs, du_inf_p, offs);
            }
            // constribution of inequality - slack constraints
            for (int k = 0; k < K; k++)
            {
                const int nu = nu_p[k];
                const int nx = nx_p[k];
                const int ng_ineq = ng_ineq_p[k];
                const int offs_g_ineq_k = offs_g_ineq[k];
                const int offs = offs_ux[k];
                GEMV_N(nu + nx, ng_ineq, 1.0, Ggt_ineq_p + k, 0, 0, lam_p, offs_g_ineq_k, 1.0, du_inf_p, offs, du_inf_p, offs);
            }
            return 0;
        }
    };
}
#endif //  DUINFEVALINCLUDED