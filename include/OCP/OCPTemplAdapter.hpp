#ifndef OCPEVALUATORINCLUDED
#define OCPEVALUATORINCLUDED
#include "OCPKKT.hpp"
#include "OCPTemplate.hpp"
#include "AUX/SmartPtr.hpp"
#include "TEMPLATES/NLP.hpp"
#define OCPMACRO1(type, name, suffix) type name##suffix = ((type)ocpkktmemory.name)
#define AUXMACRO1(type, name, suffix) type name##suffix = ((type)ocpkktmemory.aux.name)
#define SOLVERMACRO1(type, name, suffix) type name##suffix = ((type)name)
namespace fatrop
{
    class OCPTemplateAdapter : public NLP, public RefCountedObj // public OCP -> also include KKTmemory, OCPDims, ...
    {
        public:
        OCPTemplateAdapter(const RefCountPtr<OCPTemplate> &ocptempl_, MemoryAllocator &fma) : nuexpr(RefCountPtr<OCPTemplate>(ocptempl_)),nxexpr(RefCountPtr<OCPTemplate>(ocptempl_)),  ngexpr(RefCountPtr<OCPTemplate>(ocptempl_)), ocptempl(ocptempl_), ocpkktmemory(OCPDims(ocptempl_->get_horizon_length(),nuexpr, nxexpr, ngexpr), fma)
        {
        }
        int EvalHess(
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &scales_primal_vars,
            const FatropVecBF &lam,
            const FatropVecBF &scales_lam)
        {
            // horizon length
            int K = ocpkktmemory.K;
            // offsets
            int *offs_ux = (int *)ocpkktmemory.aux.ux_offs;
            int *offs_g = (int *)ocpkktmemory.aux.g_offs;
            int *offs_dyn_eq = (int *)ocpkktmemory.aux.dyn_eq_offs;
            OCPMACRO1(MAT *, RSQrqt, _p);
            OCPMACRO1(int *, nu, _p);
            SOLVERMACRO1(VEC *, primal_vars, _p);
            SOLVERMACRO1(VEC *, scales_primal_vars, _p);
            SOLVERMACRO1(VEC *, lam, _p);
            SOLVERMACRO1(VEC *, scales_lam, _p);
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
        int EvalJac(
            const FatropVecBF &primal_vars,
            const FatropVecBF &scales_primal_vars,
            const FatropVecBF &scales_lam)
        {
            // horizon length
            int K = ocpkktmemory.K;
            // offsets
            int *offs_ux = (int *)ocpkktmemory.aux.ux_offs;
            int *offs_g = (int *)ocpkktmemory.aux.g_offs;
            int *offs_dyn_eq = (int *)ocpkktmemory.aux.dyn_eq_offs;
            OCPMACRO1(MAT *, BAbt, _p);
            OCPMACRO1(MAT *, Ggt, _p);
            OCPMACRO1(int *, nu, _p);
            OCPMACRO1(int *, ng, _p);
            SOLVERMACRO1(VEC *, primal_vars, _p);
            SOLVERMACRO1(VEC *, scales_primal_vars, _p);
            SOLVERMACRO1(VEC *, scales_lam, _p);
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
        class nxExpr : public VecExpr<nxExpr, int>
        {
        public:
            nxExpr(const RefCountPtr<OCPTemplate> &parent) : parent(parent){};
            int getEl(const int ai) const { return parent->get_nxk(ai); };
            int size() const { return parent->get_horizon_length(); };

        private:
            const RefCountPtr<OCPTemplate> parent;
        };
        class nuExpr : public VecExpr<nxExpr, int>
        {
        public:
            nuExpr(const RefCountPtr<OCPTemplate> &parent) : parent(parent){};
            int getEl(const int ai) const { return parent->get_nuk(ai); };
            int size() const { return parent->get_horizon_length(); };

        private:
            const RefCountPtr<OCPTemplate> parent;
        };
        class ngExpr : public VecExpr<nxExpr, int>
        {
        public:
            ngExpr(const RefCountPtr<OCPTemplate> &parent) : parent(parent){};
            int getEl(const int ai) const { return parent->get_ngk(ai); };
            int size() const { return parent->get_horizon_length(); };

        private:
            const RefCountPtr<OCPTemplate> parent;
        };

    public:
        nxExpr nuexpr;
        nxExpr nxexpr;
        nxExpr ngexpr;
    private:
        RefCountPtr<OCPTemplate> ocptempl;
        OCPKKTMemory ocpkktmemory;
    };
} // namespace fatrop

#endif // OCPEVALUATORINCLUDED