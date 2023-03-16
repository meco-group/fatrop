#include "ocp/OCPInitializer.hpp"
using namespace fatrop;

int OCPInitializer::AdaptKKTInitial(
    OCPKKTMemory *OCP,
    const FatropVecBF &grad)
{
    // horizon length
    int K = OCP->K;
    // offsets
    int *offs_ux = (int *)OCP->aux.ux_offs.data();
    OCPMACRO(MAT *, BAbt, _p);
    OCPMACRO(MAT *, Ggt, _p);
    OCPMACRO(MAT *, RSQrqt, _p);
    OCPMACRO(int *, nu, _p);
    OCPMACRO(int *, nx, _p);
    OCPMACRO(int *, ng, _p);
    SOLVERMACRO(VEC *, grad, _p);
    // SOLVERMACRO(VEC *, s, _p);
    for (int k = 0; k < K; k++)
    {
        int nu_k = nu_p[k];
        int nx_k = nx_p[k];
        int offs_ux_k = offs_ux[k];
        fatrop_identity(nu_k + nx_k, RSQrqt_p + k, 0, 0);
        ROWIN(nu_k + nx_k, 1.0, grad_p, offs_ux_k, RSQrqt_p + k, nu_k + nx_k, 0);
    }

    for (int k = 0; k < K - 1; k++)
    {
        int nu_k = nu_p[k];
        int nx_k = nx_p[k];
        int nx_kp1 = nx_p[k + 1];
        GESE(1, nx_kp1, 0.0, BAbt_p + k, nu_k + nx_k, 0);
    }
    for (int k = 0; k < K; k++)
    {
        int nu_k = nu_p[k];
        int nx_k = nx_p[k];
        int ng_k = ng_p[k];
        if (ng_k > 0)
        {
            GESE(1, ng_k, 0.0, Ggt_p + k, nu_k + nx_k, 0);
        }
    }
    OCPMACRO(MAT *, Ggt_ineq, _p);
    OCPMACRO(int *, ng_ineq, _p);
    for (int k = 0; k < K; k++)
    {
        int nu_k = nu_p[k];
        int nx_k = nx_p[k];
        int ng_ineq_k = ng_ineq_p[k];
        if (ng_ineq_k > 0)
        {
            GESE(1, ng_ineq_k, 0.0, Ggt_ineq_p + k, nu_k + nx_k, 0);
        }
    }
    return 0;
}
int OCPInitializer::CalcSlacks(
    OCPKKTMemory *OCP,
    FatropVecBF &s)
{
    // horizon length
    int K = OCP->K;
    // offsets
    int *offs_ineq_p = (int *)OCP->aux.ineq_offs.data();
    OCPMACRO(MAT *, Ggt_ineq, _p);
    OCPMACRO(int *, nu, _p);
    OCPMACRO(int *, nx, _p);
    OCPMACRO(int *, ng_ineq, _p);
    SOLVERMACRO(VEC *, s, _p);
    for (int k = 0; k < K; k++)
    {
        int nu_k = nu_p[k];
        int nx_k = nx_p[k];
        int ng_ineq_k = ng_ineq_p[k];
        int offs_ineq_k = offs_ineq_p[k];
        if (ng_ineq_k > 0)
        {
            ROWEX(ng_ineq_k, 1.0, Ggt_ineq_p + k, nu_k + nx_k, 0, s_p, offs_ineq_k);
            // GESE(1, ng_ineq_k, 0.0, Ggt_ineq_p + k, nu_k + nx_k, 0);
        }
    }
    return 0;
}