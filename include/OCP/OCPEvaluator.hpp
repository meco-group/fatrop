#ifndef OCPEVALUATORINCLUDED
#define OCPEVALUATORINCLUDED
#include "FatropOCPKKT.hpp"
#include "FUNCTION_EVALUATION/FunctionEvaluation.hpp"
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP->name)
#define AUXMACRO(type, name, suffix) type name##suffix = ((type)OCP->aux.name)
#define SOLVERMACRO(type, name, suffix) type name##suffix = ((type)name)
namespace fatrop
{
    class OCP_evaluator
    {
        int eval_hess(OCP_KKT *OCP, double obj_scale, const fatrop_vector_bf &primal_vars, const fatrop_vector_bf &scales_primal_vars, const fatrop_vector_bf &lam, const fatrop_vector_bf &scales_lam)
        {
            // horizon length
            int K = OCP->K;
            // offsets
            int *offs_ux = (int *)OCP->aux.ux_offs;
            int *offs_g = (int *)OCP->aux.g_offs;
            int *offs_dyn_eq = (int *)OCP->aux.dyn_eq_offs;
            OCPMACRO(MAT *, RSQrqt, _p);
            OCPMACRO(int *, nu, _p);
            OCPMACRO(int *, nx, _p);
            OCPMACRO(int *, ng, _p);
            SOLVERMACRO(VEC *, primal_vars, _p);
            SOLVERMACRO(VEC *, scales_primal_vars, _p);
            SOLVERMACRO(VEC *, lam, _p);
            SOLVERMACRO(VEC *, scales_lam, _p);
            double *primal_data = primal_vars_p->pa;
            double *scales_ux_data = scales_primal_vars_p->pa;
            double *lam_data = lam_p->pa;
            double *scales_lam_data = scales_lam_p->pa;
            for (int k = 0; k < K; k++)
            {
                int nu_k = nu_p[k];
                int nu_kp1 = nu_p[k + 1];
                int offs_ux_k = offs_ux[k];
                int offs_ux_kp1 = offs_ux[k + 1];
                const double *arg[3];
                RSQrqtf.at(k)->eval_bf(arg, RSQrqt_p + k);
            }
            return 0;
        }
        int eval_jac(OCP_KKT *OCP, const fatrop_vector_bf &primal_vars, const fatrop_vector_bf &scales_primal_vars)
        {
            // horizon length
            int K = OCP->K;
            // offsets
            int *offs_ux = (int *)OCP->aux.ux_offs;
            OCPMACRO(MAT *, BAbt, _p);
            OCPMACRO(int *, nu, _p);
            SOLVERMACRO(VEC *, primal_vars, _p);
            SOLVERMACRO(VEC *, scales_primal_vars, _p);
            double *primal_data = primal_vars_p->pa;
            double *scales_primal_data = scales_primal_vars_p->pa;

            for (int k = 0; k < K - 1; k++)
            {
                int nu_k = nu_p[k];
                int nu_kp1 = nu_p[k + 1];
                int offs_ux_k = offs_ux[k];
                int offs_ux_kp1 = offs_ux[k + 1];
                const double *arg[6];
                arg[0] = primal_data + offs_ux_kp1 + nu_kp1;        // states_kp1
                arg[1] = scales_primal_data + offs_ux_kp1 + nu_kp1; // scales states_kp1
                arg[2] = primal_data + offs_ux_k + nu_k;            // states_k
                arg[3] = scales_primal_data + offs_ux_k + nu_k;     // scales states_k
                arg[4] = primal_data + offs_ux_k;                   // inputs_k
                arg[5] = scales_primal_data + offs_ux_k;            // scales inputs_k
                BAbtf.at(k)->eval_bf(arg, BAbt_p + k);
            }
            // OCPMACRO(MAT *, Ggt, _p);
            return 0;
        }

    public:
        vector<fatrop_eval_base *> BAbtf;
        vector<fatrop_eval_base *> RSQrqtf;
        vector<fatrop_eval_base *> Ggtf;
    };
} // namespace fatrop

#endif // OCPEVALUATORINCLUDED