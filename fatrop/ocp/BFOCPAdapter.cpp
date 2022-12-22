#include "ocp/BFOCPAdapter.hpp"
using namespace fatrop;
int BFOCPAdapter::evalHess(
    OCPKKTMemory *OCP,
    double obj_scale,
    const FatropVecBF &primal_vars,
    const FatropVecBF &lam)
{
    // horizon length
    int K = OCP->K;
    // offsets
    const int *offs_ux = (const int *)OCP->aux.ux_offs.data();
    int *offs_g = (int *)OCP->aux.g_offs.data();
    int *offs_dyn_eq = (int *)OCP->aux.dyn_eq_offs.data();
    int *offs_ineq = (int *)OCP->aux.g_ineq_offs.data();
    int *offs_stageparams_p = (int *)offs_stageparams.data();
    double *stageparams_p = (double *)stageparams.data();
    double *globalparams_p = (double *)globalparams.data();
    OCPMACRO(MAT *, RSQrqt, _p);
    OCPMACRO(int *, nu, _p);
    OCPMACRO(int *, nx, _p);
    SOLVERMACRO(VEC *, primal_vars, _p);
    SOLVERMACRO(VEC *, lam, _p);
    double *primal_data = primal_vars_p->pa;
    double *lam_data = lam_p->pa;

    #ifdef ENABLE_MULTITHREADING
    #pragma omp parallel for
    #endif
    for (int k = 0; k < K; k++)
    {
        int nu_k = nu_p[k];
        int nx_k = nx_p[k];
        int offs_ux_k = offs_ux[k];
        int offs_dyn_eq_k = offs_dyn_eq[k];
        int offs_g_k = offs_g[k];
        int offs_ineq_k = offs_ineq[k];
        int offs_stageparams_k = offs_stageparams_p[k];
        ocptempl->eval_RSQrqtk(
            &obj_scale,
            primal_data + offs_ux_k,
            primal_data + offs_ux_k + nu_k,
            lam_data + offs_dyn_eq_k,
            lam_data + offs_g_k,
            lam_data + offs_ineq_k,
            stageparams_p + offs_stageparams_k,
            globalparams_p,
            RSQrqt_p + k,
            k);
        if (k > 0)
        {
            int offs_dyn_eq_km1 = offs_dyn_eq[k - 1];
            ROWAD(nx_k, -1.0, lam_p, offs_dyn_eq_km1, RSQrqt_p + k, nu_k + nx_k, nu_k);
        }
        // std::cout << "k " << k << std::endl;
        // blasfeo_print_dmat(nu_k+nx_k+1, nu_k+nx_k,  RSQrqt_p +k, 0,0);
    }
    return 0;
}
int BFOCPAdapter::evalJac(
    OCPKKTMemory *OCP,
    const FatropVecBF &primal_vars,
    const FatropVecBF &slack_vars)
{
    // horizon length
    int K = OCP->K;
    // offsets
    int *offs_ux = (int *)OCP->aux.ux_offs.data();
    int *offs_ineq = (int *)OCP->aux.ineq_offs.data();
    int *offs_stageparams_p = (int *)offs_stageparams.data();
    double *stageparams_p = (double *)stageparams.data();
    double *globalparams_p = (double *)globalparams.data();
    OCPMACRO(MAT *, BAbt, _p);
    OCPMACRO(MAT *, Ggt, _p);
    OCPMACRO(MAT *, Ggt_ineq, _p);
    OCPMACRO(int *, nu, _p);
    OCPMACRO(int *, nx, _p);
    OCPMACRO(int *, ng, _p);
    OCPMACRO(int *, ng_ineq, _p);
    SOLVERMACRO(VEC *, primal_vars, _p);
    double *primal_data = primal_vars_p->pa;

    #ifdef ENABLE_MULTITHREADING
    #pragma omp parallel for
    #endif
    for (int k = 0; k < K - 1; k++)
    {
        int nu_k = nu_p[k];
        int nu_kp1 = nu_p[k + 1];
        int offs_ux_k = offs_ux[k];
        int offs_ux_kp1 = offs_ux[k + 1];
        int offs_stageparams_k = offs_stageparams_p[k];
        ocptempl->eval_BAbtk(
            primal_data + offs_ux_kp1 + nu_kp1,
            primal_data + offs_ux_k,
            primal_data + offs_ux_k + nu_k,
            stageparams_p + offs_stageparams_k,
            globalparams_p,
            BAbt_p + k,
            k);
    }
    #ifdef ENABLE_MULTITHREADING
    #pragma omp parallel for
    #endif
    for (int k = 0; k < K; k++)
    {
        int nu_k = nu_p[k];
        int ng_k = ng_p[k];
        int offs_ux_k = offs_ux[k];
        int offs_stageparams_k = offs_stageparams_p[k];
        if (ng_k > 0)
        {
            ocptempl->eval_Ggtk(
                primal_data + offs_ux_k,
                primal_data + offs_ux_k + nu_k,
                stageparams_p + offs_stageparams_k,
                globalparams_p,
                Ggt_p + k,
                k);
        }
    }
    VEC *slack_vars_bf = (VEC *)slack_vars;
    #ifdef ENABLE_MULTITHREADING
    #pragma omp parallel for
    #endif
    for (int k = 0; k < K; k++)
    {
        int nu_k = nu_p[k];
        int nx_k = nx_p[k];
        int ng_ineq_k = ng_ineq_p[k];
        int offs_ux_k = offs_ux[k];
        int offs_ineq_k = offs_ineq[k];
        int offs_stageparams_k = offs_stageparams_p[k];
        if (ng_ineq_k > 0)
        {
            ocptempl->eval_Ggt_ineqk(
                primal_data + offs_ux_k,
                primal_data + offs_ux_k + nu_k,
                stageparams_p + offs_stageparams_k,
                globalparams_p,
                Ggt_ineq_p + k,
                k);
            // rewrite problem
            ROWAD(ng_ineq_k, -1.0, slack_vars_bf, offs_ineq_k, Ggt_ineq_p + k, nu_k + nx_k, 0);
        }
    }
    return 0;
}
int BFOCPAdapter::EvalConstraintViolation(
    OCPKKTMemory *OCP,
    const FatropVecBF &primal_vars,
    const FatropVecBF &slack_vars,
    FatropVecBF &constraint_violation)
{
    // horizon length
    int K = OCP->K;
    // offsets
    int *offs_ux = (int *)OCP->aux.ux_offs.data();
    int *offs_g = (int *)OCP->aux.g_offs.data();
    int *offs_dyn_eq = (int *)OCP->aux.dyn_eq_offs.data();
    int *offs_ineq = (int *)OCP->aux.ineq_offs.data();
    int *offs_g_ineq = (int *)OCP->aux.g_ineq_offs.data();
    int *offs_stageparams_p = (int *)offs_stageparams.data();
    double *stageparams_p = (double *)stageparams.data();
    double *globalparams_p = (double *)globalparams.data();
    double *cv_p = ((VEC *)constraint_violation)->pa;
    OCPMACRO(int *, nu, _p);
    OCPMACRO(int *, ng, _p);
    OCPMACRO(int *, ng_ineq, _p);
    SOLVERMACRO(VEC *, primal_vars, _p);
    double *primal_data = primal_vars_p->pa;
#ifdef ENABLE_MULTITHREADING
#pragma omp parallel for
#endif
    for (int k = 0; k < K - 1; k++)
    {
        int nu_k = nu_p[k];
        int nu_kp1 = nu_p[k + 1];
        int offs_ux_k = offs_ux[k];
        int offs_ux_kp1 = offs_ux[k + 1];
        int offs_dyn_eq_k = offs_dyn_eq[k];
        int offs_stageparams_k = offs_stageparams_p[k];
        ocptempl->eval_bk(
            primal_data + offs_ux_kp1 + nu_kp1,
            primal_data + offs_ux_k,
            primal_data + offs_ux_k + nu_k,
            stageparams_p + offs_stageparams_k,
            globalparams_p,
            cv_p + offs_dyn_eq_k,
            k);
    }
    #ifdef ENABLE_MULTITHREADING
    #pragma omp parallel for
    #endif
    for (int k = 0; k < K; k++)
    {
        int ng_k = ng_p[k];
        if (ng_k > 0)
        {
            int nu_k = nu_p[k];
            int offs_ux_k = offs_ux[k];
            int offs_g_k = offs_g[k];
            int offs_stageparams_k = offs_stageparams_p[k];
            ocptempl->eval_gk(
                primal_data + offs_ux_k,
                primal_data + offs_ux_k + nu_k,
                stageparams_p + offs_stageparams_k,
                globalparams_p,
                cv_p + offs_g_k,
                k);
        }
    }
    VEC *cv_bf = (VEC *)constraint_violation;
    VEC *slack_vars_bf = (VEC *)slack_vars;
    #ifdef ENABLE_MULTITHREADING
    #pragma omp parallel for
    #endif
    for (int k = 0; k < K; k++)
    {
        int ng_ineq_k = ng_ineq_p[k];
        if (ng_ineq_k > 0)
        {
            int nu_k = nu_p[k];
            int ng_ineq_k = ng_ineq_p[k];
            int offs_ux_k = offs_ux[k];
            int offs_ineq_k = offs_ineq[k];
            int offs_g_ineq_k = offs_g_ineq[k];
            int offs_stageparams_k = offs_stageparams_p[k];
            ocptempl->eval_gineqk(
                primal_data + offs_ux_k,
                primal_data + offs_ux_k + nu_k,
                stageparams_p + offs_stageparams_k,
                globalparams_p,
                cv_p + offs_g_ineq_k,
                k);
            // rewrite problem
            AXPY(ng_ineq_k, -1.0, slack_vars_bf, offs_ineq_k, cv_bf, offs_g_ineq_k, cv_bf, offs_g_ineq_k);
        }
    }
    return 0;
}
int BFOCPAdapter::EvalGrad(
    OCPKKTMemory *OCP,
    double obj_scale,
    const FatropVecBF &primal_vars,
    FatropVecBF &gradient)
{
    // horizon length
    int K = OCP->K;
    // offsets
    const int *offs_ux = (const int *)OCP->aux.ux_offs.data();
    double *grad_p = ((VEC *)gradient)->pa;
    OCPMACRO(int *, nu, _p);
    SOLVERMACRO(VEC *, primal_vars, _p);
    double *primal_data = primal_vars_p->pa;
    int *offs_stageparams_p = (int *)offs_stageparams.data();
    double *stageparams_p = (double *)stageparams.data();
    double *globalparams_p = (double *)globalparams.data();
    #ifdef ENABLE_MULTITHREADING
    #pragma omp parallel for
    #endif
    for (int k = 0; k < K; k++)
    {
        int nu_k = nu_p[k];
        int offs_ux_k = offs_ux[k];
        int offs_stageparams_k = offs_stageparams_p[k];
        ocptempl->eval_rqk(
            &obj_scale,
            primal_data + offs_ux_k,
            primal_data + offs_ux_k + nu_k,
            stageparams_p + offs_stageparams_k,
            globalparams_p,
            grad_p + offs_ux_k,
            k);
    }
    return 0;
};
int BFOCPAdapter::EvalObj(
    OCPKKTMemory *OCP,
    double obj_scale,
    const FatropVecBF &primal_vars,
    double &res)
{
    // horizon length
    int K = OCP->K;
    // offsets
    const int *offs_ux = (const int *)OCP->aux.ux_offs.data();
    int *offs_stageparams_p = (int *)offs_stageparams.data();
    double *stageparams_p = (double *)stageparams.data();
    double *globalparams_p = (double *)globalparams.data();
    OCPMACRO(int *, nu, _p);
    SOLVERMACRO(VEC *, primal_vars, _p);
    double *primal_data = primal_vars_p->pa;
    double restot = 0.0;
    for (int k = 0; k < K; k++)
    {
        int nu_k = nu_p[k];
        int offs_ux_k = offs_ux[k];
        int offs_stageparams_k = offs_stageparams_p[k];
        double resk = 0.0;
        ocptempl->eval_Lk(
            &obj_scale,
            primal_data + offs_ux_k,
            primal_data + offs_ux_k + nu_k,
            stageparams_p + offs_stageparams_k,
            globalparams_p,
            &resk,
            k);
        restot += resk;
    }
    res = restot;
    return 0;
};

int BFOCPAdapter::EvalDynamics(
    OCPKKTMemory *OCP,
    const int k,
    const FatropVecBF &uk,
    const FatropVecBF &xk,
    FatropVecBF &xkp1)
{
    // offsets
    int *offs_stageparams_p = (int *)offs_stageparams.data();
    int offs_stageparams_k = offs_stageparams_p[k];
    double *stageparams_p = (double *)stageparams.data();
    double *globalparams_p = (double *)globalparams.data();
    double *ukp = ((blasfeo_dvec *)uk)->pa + uk.offset();
    double *xkp = ((blasfeo_dvec *)xk)->pa + xk.offset();
    double *xkp1p = ((blasfeo_dvec *)xkp1)->pa + xkp1.offset();
    double *x_dummy_p = x_dummy.data();
    ocptempl->eval_bk(
        x_dummy_p,
        ukp,
        xkp,
        stageparams_p + offs_stageparams_k,
        globalparams_p,
        xkp1p,
        k);
    return 0;
};
void BFOCPAdapter::SetParams(const vector<double> &stage_params_in, const vector<double> &global_params_in)
{
    stageparams = stage_params_in;
    globalparams = global_params_in;
    return;
}
void BFOCPAdapter::SetInitial(const shared_ptr<FatropData> &fatropdata, vector<double> &initial_u, vector<double> &initial_x)
{
    // offsets
    VEC *ux_intial_p = (VEC *)fatropdata->x_initial;
    double *u_p = initial_u.data();
    double *x_p = initial_x.data();
    int offs_nu = 0;
    int offs_nx = 0;
    int offs_nux = 0;
    for (int k = 0; k < K; k++)
    {
        int nu_k = ocptempl->get_nuk(k);
        int nx_k = ocptempl->get_nxk(k);
        PACKVEC(nu_k, u_p+ offs_nu, 1, ux_intial_p, offs_nux);
        PACKVEC(nx_k, x_p+ offs_nx, 1, ux_intial_p, offs_nux + nu_k);
        offs_nu += nu_k;
        offs_nx += nx_k;
        offs_nux += nu_k + nx_k;
    }
    return;
}
void BFOCPAdapter::GetSolution(const shared_ptr<FatropData> &fatropdata, vector<double> &u, vector<double> &x)
{
    // offsets
    VEC *ux_sol = (VEC *)fatropdata->x_curr;
    double *u_p = u.data();
    double *x_p = x.data();
    int offs_nu = 0;
    int offs_nx = 0;
    int offs_nux = 0;
    for (int k = 0; k < K; k++)
    {
        int nu_k = ocptempl->get_nuk(k);
        int nx_k = ocptempl->get_nxk(k);
        UNPACKVEC(nu_k, ux_sol, offs_nux, u_p+ offs_nu, 1);
        UNPACKVEC(nx_k, ux_sol, offs_nux + nu_k, x_p+ offs_nx, 1);
        offs_nu += nu_k;
        offs_nx += nx_k;
        offs_nux += nu_k + nx_k;
    }
    return;
}