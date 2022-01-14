#ifndef OCPEVALUATORINCLUDED
#define OCPEVALUATORINCLUDED
#include "OCPKKT.hpp"
#include "OCPTemplate.hpp"
#include "AUX/SmartPtr.hpp"
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP->name)
#define AUXMACRO(type, name, suffix) type name##suffix = ((type)OCP->aux.name)
#define SOLVERMACRO(type, name, suffix) type name##suffix = ((type)name)
namespace fatrop
{
    class OCPTemplateAdapter // public OCP -> also include KKTmemory, OCPDims, ...
    {
        OCPTemplateAdapter(const RefCountPtr<OCPTemplate>& ocptempl):ocptempl(ocptempl)
        {
        }
        int eval_hess(OCPKKTMemory *OCP,
                      double obj_scale,
                      const FatropVecBF &primal_vars,
                      const FatropVecBF &scales_primal_vars,
                      const FatropVecBF &lam,
                      const FatropVecBF &scales_lam)
        {
            // horizon length
            int K = OCP->K;
            // offsets
            int *offs_ux = (int *)OCP->aux.ux_offs;
            int *offs_g = (int *)OCP->aux.g_offs;
            int *offs_dyn_eq = (int *)OCP->aux.dyn_eq_offs;
            OCPMACRO(MAT *, RSQrqt, _p);
            OCPMACRO(int *, nu, _p);
            SOLVERMACRO(VEC *, primal_vars, _p);
            SOLVERMACRO(VEC *, scales_primal_vars, _p);
            SOLVERMACRO(VEC *, lam, _p);
            SOLVERMACRO(VEC *, scales_lam, _p);
            double *primal_data = primal_vars_p->pa;
            double *scales_primal_data = scales_primal_vars_p->pa;
            double *lam_data = lam_p->pa;
            double *scales_lam_data = scales_lam_p->pa;
            for (int k = 0; k < K; k++)
            {
                int nu_k = nu_p[k];
                int offs_ux_k = offs_ux[k];
                int offs_dyn_eq_k = offs_dyn_eq[k];
                int offs_g_k = offs_g[k];
                ocptempl->eval_RSQrqtk(
                    &obj_scale,
                    primal_data + offs_ux_k,
                    scales_primal_data + offs_ux_k,
                    primal_data + offs_ux_k + nu_k,
                    scales_primal_data + offs_ux_k + nu_k,
                    lam_data + offs_dyn_eq_k,
                    scales_lam_data + offs_dyn_eq_k,
                    lam_data + offs_g_k,
                    scales_lam_data + offs_g_k,
                    RSQrqt_p + k,
                    k);
            }
            return 0;
        }
        int eval_jac(OCPKKTMemory *OCP,
                     const FatropVecBF &primal_vars,
                     const FatropVecBF &scales_primal_vars,
                     const FatropVecBF &scales_lam)
        {
            // horizon length
            int K = OCP->K;
            // offsets
            int *offs_ux = (int *)OCP->aux.ux_offs;
            int *offs_g = (int *)OCP->aux.g_offs;
            int *offs_dyn_eq = (int *)OCP->aux.dyn_eq_offs;
            OCPMACRO(MAT *, BAbt, _p);
            OCPMACRO(MAT *, Ggt, _p);
            OCPMACRO(int *, nu, _p);
            OCPMACRO(int *, ng, _p);
            SOLVERMACRO(VEC *, primal_vars, _p);
            SOLVERMACRO(VEC *, scales_primal_vars, _p);
            SOLVERMACRO(VEC *, scales_lam, _p);
            double *primal_data = primal_vars_p->pa;
            double *scales_primal_data = scales_primal_vars_p->pa;
            double *scales_lam_data = scales_lam_p->pa;

            for (int k = 0; k < K - 1; k++)
            {
                int nu_k = nu_p[k];
                int ng_k = ng_p[k];
                int nu_kp1 = nu_p[k + 1];
                int offs_ux_k = offs_ux[k];
                int offs_ux_kp1 = offs_ux[k + 1];
                int offs_dyn_eq_k = offs_dyn_eq[k];
                int offs_g_k = offs_g[k];
                ocptempl->eval_BAbtk(
                    primal_data + offs_ux_kp1 + nu_kp1,
                    scales_primal_data + offs_ux_kp1 + nu_kp1,
                    primal_data + offs_ux_k + nu_k,
                    scales_primal_data + offs_ux_k + nu_k,
                    primal_data + offs_ux_k,
                    scales_primal_data + offs_ux_k,
                    scales_lam_data + offs_dyn_eq_k,
                    BAbt_p + k,
                    k);
                if (ng_k > 0)
                {
                    ocptempl->eval_Ggtk(
                        primal_data + offs_ux_k + nu_k,
                        scales_primal_data + offs_ux_k + nu_k,
                        primal_data + offs_ux_k,
                        scales_primal_data + offs_ux_k,
                        scales_lam_data + offs_g_k,
                        Ggt_p + k,
                        k);
                }
            }
            return 0;
        }
        int eval_g()
        {
            return 0;
        }

    private:
        RefCountPtr<OCPTemplate> ocptempl;
    };
} // namespace fatrop

#endif // OCPEVALUATORINCLUDED