#ifndef OCPEVALUATORINCLUDED
#define OCPEVALUATORINCLUDED
#include "OCPKKT.hpp"
#include "BFOCP.hpp"
#include "AUX/SmartPtr.hpp"
#include "OCP.hpp"
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP->name)
#define AUXMACRO(type, name, suffix) type name##suffix = ((type)OCP->aux.name)
#define SOLVERMACRO(type, name, suffix) type name##suffix = ((type)name)
namespace fatrop
{
    class BFOCPAdapter : public OCP // public OCP -> also include KKTmemory, OCPDims, ...
    {
    public:
        BFOCPAdapter(const RefCountPtr<BFOCP> &ocptempl_) : nuexpr(RefCountPtr<BFOCP>(ocptempl_)), nxexpr(RefCountPtr<BFOCP>(ocptempl_)), ngexpr(RefCountPtr<BFOCP>(ocptempl_)), ocptempl(ocptempl_)
        {
        }
        int evalHess(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &lam) override
        {
            // horizon length
            int K = OCP->K;
            // offsets
            const int *offs_ux = (const int *)OCP->aux.ux_offs.data();
            int *offs_g = (int *)OCP->aux.g_offs.data();
            int *offs_dyn_eq = (int *)OCP->aux.dyn_eq_offs.data();
            OCPMACRO(MAT *, RSQrqt, _p);
            OCPMACRO(int *, nu, _p);
            SOLVERMACRO(VEC *, primal_vars, _p);
            SOLVERMACRO(VEC *, lam, _p);
            double *primal_data = primal_vars_p->pa;
            double *lam_data = lam_p->pa;
            for (int k = 0; k < K; k++)
            {
                int nu_k = nu_p[k];
                int offs_ux_k = offs_ux[k];
                int offs_dyn_eq_k = offs_dyn_eq[k];
                int offs_g_k = offs_g[k];
                ocptempl->eval_RSQrqtk(
                    &obj_scale,
                    primal_data + offs_ux_k,
                    primal_data + offs_ux_k + nu_k,
                    lam_data + offs_dyn_eq_k,
                    lam_data + offs_g_k,
                    RSQrqt_p + k,
                    k);
            }
            return 0;
        }
        int evalJac(
            OCPKKTMemory *OCP,
            const FatropVecBF &primal_vars ) override
        {
            // horizon length
            int K = OCP->K;
            // offsets
            int *offs_ux = (int *)OCP->aux.ux_offs.data();
            OCPMACRO(MAT *, BAbt, _p);
            OCPMACRO(MAT *, Ggt, _p);
            OCPMACRO(int *, nu, _p);
            OCPMACRO(int *, ng, _p);
            SOLVERMACRO(VEC *, primal_vars, _p);
            double *primal_data = primal_vars_p->pa;

            for (int k = 0; k < K - 1; k++)
            {
                int nu_k = nu_p[k];
                int ng_k = ng_p[k];
                int nu_kp1 = nu_p[k + 1];
                int offs_ux_k = offs_ux[k];
                int offs_ux_kp1 = offs_ux[k + 1];
                ocptempl->eval_BAbtk(
                    primal_data + offs_ux_kp1 + nu_kp1,
                    primal_data + offs_ux_k + nu_k,
                    primal_data + offs_ux_k,
                    BAbt_p + k,
                    k);
                if (ng_k > 0)
                {
                    ocptempl->eval_Ggtk(
                        primal_data + offs_ux_k + nu_k,
                        primal_data + offs_ux_k,
                        Ggt_p + k,
                        k);
                }
            }
            {
                int nu_k = nu_p[K - 1];
                int ng_k = ng_p[K - 1];
                int offs_ux_k = offs_ux[K - 1];
                if (ng_k > 0)
                {
                    ocptempl->eval_Ggtk(
                        primal_data + offs_ux_k + nu_k,
                        primal_data + offs_ux_k,
                        Ggt_p + K - 1,
                        K - 1);
                }
            }
            return 0;
        }
        int EvalConstraintViolation(
            OCPKKTMemory *OCP,
            const FatropVecBF &primal_vars,
            FatropVecBF &constraint_violation) override
        {
            // horizon length
            int K = OCP->K;
            // offsets
            int *offs_ux = (int *)OCP->aux.ux_offs.data();
            int *offs_g = (int *)OCP->aux.g_offs.data();
            int *offs_dyn_eq = (int *)OCP->aux.dyn_eq_offs.data();
            double *cv_p = ((VEC *)constraint_violation)->pa;
            OCPMACRO(int *, nu, _p);
            OCPMACRO(int *, ng, _p);
            SOLVERMACRO(VEC *, primal_vars, _p);
            double *primal_data = primal_vars_p->pa;

            for (int k = 0; k < K - 1; k++)
            {
                int nu_k = nu_p[k];
                int ng_k = ng_p[k];
                int nu_kp1 = nu_p[k + 1];
                int offs_ux_k = offs_ux[k];
                int offs_ux_kp1 = offs_ux[k + 1];
                int offs_dyn_eq_k = offs_dyn_eq[k];
                int offs_g_k = offs_g[k];
                ocptempl->eval_bk(
                    primal_data + offs_ux_kp1 + nu_kp1,
                    primal_data + offs_ux_k + nu_k,
                    primal_data + offs_ux_k,
                    cv_p + offs_dyn_eq_k,
                    k);
                if (ng_k > 0)
                {
                    ocptempl->eval_gk(
                        primal_data + offs_ux_k + nu_k,
                        primal_data + offs_ux_k,
                        cv_p + offs_g_k,
                        k);
                }
            }
            {
                int nu_k = nu_p[K - 1];
                int ng_k = ng_p[K - 1];
                int offs_ux_k = offs_ux[K - 1];
                int offs_g_k = offs_g[K - 1];
                if (ng_k > 0)
                {
                    ocptempl->eval_gk(
                        primal_data + offs_ux_k + nu_k,
                        primal_data + offs_ux_k,
                        cv_p + offs_g_k,
                        K - 1);
                }
            }
            return 0;
        }
        int EvalGrad(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            FatropVecBF &gradient) override
        {
            // horizon length
            int K = OCP->K;
            // offsets
            const int *offs_ux = (const int *)OCP->aux.ux_offs.data();
            double *grad_p = ((VEC *)gradient)->pa;
            OCPMACRO(int *, nu, _p);
            SOLVERMACRO(VEC *, primal_vars, _p);
            double *primal_data = primal_vars_p->pa;
            for (int k = 0; k < K; k++)
            {
                int nu_k = nu_p[k];
                int offs_ux_k = offs_ux[k];
                ocptempl->eval_rqk(
                    &obj_scale,
                    primal_data + offs_ux_k,
                    primal_data + offs_ux_k + nu_k,
                    grad_p + offs_ux_k,
                    k);
            }
            return 0;
        };
        int EvalObj(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            double &res)
        {
            // horizon length
            int K = OCP->K;
            // offsets
            const int *offs_ux = (const int *)OCP->aux.ux_offs.data();
            OCPMACRO(int *, nu, _p);
            SOLVERMACRO(VEC *, primal_vars, _p);
            double *primal_data = primal_vars_p->pa;
            double restot = 0.0;
            for (int k = 0; k < K; k++)
            {
                int nu_k = nu_p[k];
                int offs_ux_k = offs_ux[k];
                double resk = 0.0;
                ocptempl->eval_rqk(
                    &obj_scale,
                    primal_data + offs_ux_k,
                    primal_data + offs_ux_k + nu_k,
                    &resk,
                    k);
                restot += resk;
            }
            res = restot;
            return 0;
        }

        ;
        OCPDims GetOCPDims() const override
        {
            return OCPDims(ocptempl->get_horizon_length(), nuexpr, nxexpr, ngexpr);
        }

    private:
        class nxExpr : public VecExpr<nxExpr, int>
        {
        public:
            nxExpr(const RefCountPtr<BFOCP> &parent) : parent(parent){};
            int getEl(const int ai) const { return parent->get_nxk(ai); };
            int size() const { return parent->get_horizon_length(); };

        private:
            const RefCountPtr<BFOCP> parent;
        };
        class nuExpr : public VecExpr<nuExpr, int>
        {
        public:
            nuExpr(const RefCountPtr<BFOCP> &parent) : parent(parent){};
            int getEl(const int ai) const { return parent->get_nuk(ai); };
            int size() const { return parent->get_horizon_length(); };

        private:
            const RefCountPtr<BFOCP> parent;
        };
        class ngExpr : public VecExpr<ngExpr, int>
        {
        public:
            ngExpr(const RefCountPtr<BFOCP> &parent) : parent(parent){};
            int getEl(const int ai) const { return parent->get_ngk(ai); };
            int size() const { return parent->get_horizon_length(); };

        private:
            const RefCountPtr<BFOCP> parent;
        };

    public:
        nuExpr nuexpr;
        nxExpr nxexpr;
        ngExpr ngexpr;

    private:
        RefCountPtr<BFOCP> ocptempl;
    };
} // namespace fatrop

#endif // OCPEVALUATORINCLUDED