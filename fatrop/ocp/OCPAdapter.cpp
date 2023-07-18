/*
 * Fatrop - A fast trajectory optimization solver
 * Copyright (C) 2022, 2023 Lander Vanroye, KU Leuven. All rights reserved.
 *
 * This file is part of Fatrop.
 *
 * Fatrop is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fatrop is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Fatrop.  If not, see <http://www.gnu.org/licenses/>. */
#include "ocp/OCPAdapter.hpp"

using namespace fatrop;
using namespace std;
fatrop_int OCPAdapter::eval_lag_hess(
    OCPKKTMemory *OCP,
    double obj_scale,
    const FatropVecBF &primal_vars,
    const FatropVecBF &lam)
{
    // horizon length
    fatrop_int K = OCP->K;
    // offsets
    const fatrop_int *offs_ux = (const fatrop_int *)OCP->aux.ux_offs.data();
    fatrop_int *offs_g = (fatrop_int *)OCP->aux.g_offs.data();
    fatrop_int *offs_dyn_eq = (fatrop_int *)OCP->aux.dyn_eq_offs.data();
    fatrop_int *offs_ineq = (fatrop_int *)OCP->aux.g_ineq_offs.data();
    fatrop_int *offs_stageparams_p = (fatrop_int *)offs_stageparams.data();
    double *stageparams_p = (double *)stageparams.data();
    double *globalparams_p = (double *)globalparams.data();
    OCPMACRO(MAT *, RSQrqt, _p);
    OCPMACRO(fatrop_int *, nu, _p);
    OCPMACRO(fatrop_int *, nx, _p);
    SOLVERMACRO(VEC *, primal_vars, _p);
    SOLVERMACRO(VEC *, lam, _p);
    double *primal_data = primal_vars_p->pa;
    double *lam_data = lam_p->pa;

    OCPMACRO(MAT *, BAbt, _p);
    OCPMACRO(MAT *, Ggt, _p);
    OCPMACRO(MAT *, Ggt_ineq, _p);

#ifdef ENABLE_MULTITHREADING
#pragma omp parallel for
#endif
    for (fatrop_int k = 0; k < K; k++)
    {
        fatrop_int nu_k = nu_p[k];
        fatrop_int nx_k = nx_p[k];
        fatrop_int offs_ux_k = offs_ux[k];
        fatrop_int offs_dyn_eq_k = offs_dyn_eq[k];
        fatrop_int offs_g_k = offs_g[k];
        fatrop_int offs_ineq_k = offs_ineq[k];
        fatrop_int offs_stageparams_k = offs_stageparams_p[k];
        int ret = 0;
        if (bfgs)
            ret = OCPBFGS_updaters[k].update(RSQrqt_p + k, primal_vars_p, offs_ux_k, gradbuf[k], 0, BAbt_p + k, lam_p, offs_dyn_eq_k, Ggt_p + k, lam_p, offs_g_k, Ggt_ineq_p + k, lam_p, offs_ineq_k);
        if (!bfgs || ret == 1)
        {
        // std::cout << "using exact Hess " << k << std::endl;
            ret = ocptempl->eval_RSQrqtk(
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
        }

        if (k > 0)
        {
            fatrop_int offs_dyn_eq_km1 = offs_dyn_eq[k - 1];
            ROWAD(nx_k, -1.0, lam_p, offs_dyn_eq_km1, RSQrqt_p + k, nu_k + nx_k, nu_k);
        }
        // std::cout << "k " << k << std::endl;
        // blasfeo_print_dmat(nu_k+nx_k+1, nu_k+nx_k,  RSQrqt_p +k, 0,0);
    }
    return 0;
}
fatrop_int OCPAdapter::eval_constr_jac(
    OCPKKTMemory *OCP,
    const FatropVecBF &primal_vars,
    const FatropVecBF &slack_vars)
{
    // horizon length
    fatrop_int K = OCP->K;
    // offsets
    fatrop_int *offs_ux = (fatrop_int *)OCP->aux.ux_offs.data();
    fatrop_int *offs_ineq = (fatrop_int *)OCP->aux.ineq_offs.data();
    fatrop_int *offs_stageparams_p = (fatrop_int *)offs_stageparams.data();
    double *stageparams_p = (double *)stageparams.data();
    double *globalparams_p = (double *)globalparams.data();
    OCPMACRO(MAT *, BAbt, _p);
    OCPMACRO(MAT *, Ggt, _p);
    OCPMACRO(MAT *, Ggt_ineq, _p);
    OCPMACRO(fatrop_int *, nu, _p);
    OCPMACRO(fatrop_int *, nx, _p);
    OCPMACRO(fatrop_int *, ng, _p);
    OCPMACRO(fatrop_int *, ng_ineq, _p);
    SOLVERMACRO(VEC *, primal_vars, _p);
    double *primal_data = primal_vars_p->pa;

#ifdef ENABLE_MULTITHREADING
#pragma omp parallel for
#endif
    for (fatrop_int k = 0; k < K - 1; k++)
    {
        fatrop_int nu_k = nu_p[k];
        fatrop_int nu_kp1 = nu_p[k + 1];
        fatrop_int offs_ux_k = offs_ux[k];
        fatrop_int offs_ux_kp1 = offs_ux[k + 1];
        fatrop_int offs_stageparams_k = offs_stageparams_p[k];
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
    for (fatrop_int k = 0; k < K; k++)
    {
        fatrop_int nu_k = nu_p[k];
        fatrop_int ng_k = ng_p[k];
        fatrop_int offs_ux_k = offs_ux[k];
        fatrop_int offs_stageparams_k = offs_stageparams_p[k];
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
    for (fatrop_int k = 0; k < K; k++)
    {
        fatrop_int nu_k = nu_p[k];
        fatrop_int nx_k = nx_p[k];
        fatrop_int ng_ineq_k = ng_ineq_p[k];
        fatrop_int offs_ux_k = offs_ux[k];
        fatrop_int offs_ineq_k = offs_ineq[k];
        fatrop_int offs_stageparams_k = offs_stageparams_p[k];
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
fatrop_int OCPAdapter::eval_contr_viol(
    OCPKKTMemory *OCP,
    const FatropVecBF &primal_vars,
    const FatropVecBF &slack_vars,
    FatropVecBF &constraint_violation)
{
    // horizon length
    fatrop_int K = OCP->K;
    // offsets
    fatrop_int *offs_ux = (fatrop_int *)OCP->aux.ux_offs.data();
    fatrop_int *offs_g = (fatrop_int *)OCP->aux.g_offs.data();
    fatrop_int *offs_dyn_eq = (fatrop_int *)OCP->aux.dyn_eq_offs.data();
    fatrop_int *offs_ineq = (fatrop_int *)OCP->aux.ineq_offs.data();
    fatrop_int *offs_g_ineq = (fatrop_int *)OCP->aux.g_ineq_offs.data();
    fatrop_int *offs_stageparams_p = (fatrop_int *)offs_stageparams.data();
    double *stageparams_p = (double *)stageparams.data();
    double *globalparams_p = (double *)globalparams.data();
    double *cv_p = ((VEC *)constraint_violation)->pa;
    OCPMACRO(fatrop_int *, nu, _p);
    OCPMACRO(fatrop_int *, ng, _p);
    OCPMACRO(fatrop_int *, ng_ineq, _p);
    SOLVERMACRO(VEC *, primal_vars, _p);
    double *primal_data = primal_vars_p->pa;
#ifdef ENABLE_MULTITHREADING
#pragma omp parallel for
#endif
    for (fatrop_int k = 0; k < K - 1; k++)
    {
        fatrop_int nu_k = nu_p[k];
        fatrop_int nu_kp1 = nu_p[k + 1];
        fatrop_int offs_ux_k = offs_ux[k];
        fatrop_int offs_ux_kp1 = offs_ux[k + 1];
        fatrop_int offs_dyn_eq_k = offs_dyn_eq[k];
        fatrop_int offs_stageparams_k = offs_stageparams_p[k];
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
    for (fatrop_int k = 0; k < K; k++)
    {
        fatrop_int ng_k = ng_p[k];
        if (ng_k > 0)
        {
            fatrop_int nu_k = nu_p[k];
            fatrop_int offs_ux_k = offs_ux[k];
            fatrop_int offs_g_k = offs_g[k];
            fatrop_int offs_stageparams_k = offs_stageparams_p[k];
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
    for (fatrop_int k = 0; k < K; k++)
    {
        fatrop_int ng_ineq_k = ng_ineq_p[k];
        if (ng_ineq_k > 0)
        {
            fatrop_int nu_k = nu_p[k];
            fatrop_int ng_ineq_k = ng_ineq_p[k];
            fatrop_int offs_ux_k = offs_ux[k];
            fatrop_int offs_ineq_k = offs_ineq[k];
            fatrop_int offs_g_ineq_k = offs_g_ineq[k];
            fatrop_int offs_stageparams_k = offs_stageparams_p[k];
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
fatrop_int OCPAdapter::eval_obj_grad(
    OCPKKTMemory *OCP,
    double obj_scale,
    const FatropVecBF &primal_vars,
    FatropVecBF &gradient)
{
    // horizon length
    fatrop_int K = OCP->K;
    // offsets
    const fatrop_int *offs_ux = (const fatrop_int *)OCP->aux.ux_offs.data();
    double *grad_p = ((VEC *)gradient)->pa;
    OCPMACRO(fatrop_int *, nu, _p);
    OCPMACRO(fatrop_int *, nx, _p);
    SOLVERMACRO(VEC *, primal_vars, _p);
    SOLVERMACRO(VEC *, gradient, _p);
    double *primal_data = primal_vars_p->pa;
    fatrop_int *offs_stageparams_p = (fatrop_int *)offs_stageparams.data();
    double *stageparams_p = (double *)stageparams.data();
    double *globalparams_p = (double *)globalparams.data();
#ifdef ENABLE_MULTITHREADING
#pragma omp parallel for
#endif
    for (fatrop_int k = 0; k < K; k++)
    {
        fatrop_int nu_k = nu_p[k];
        fatrop_int nx_k = nx_p[k];
        fatrop_int offs_ux_k = offs_ux[k];
        fatrop_int offs_stageparams_k = offs_stageparams_p[k];
        ocptempl->eval_rqk(
            &obj_scale,
            primal_data + offs_ux_k,
            primal_data + offs_ux_k + nu_k,
            stageparams_p + offs_stageparams_k,
            globalparams_p,
            grad_p + offs_ux_k,
            k);
        // save result in grad buf
        VECCP(nu_k + nx_k, gradient_p, offs_ux_k, gradbuf[k], 0);
        // blasfeo_print_dvec(nu_k + nx_k, gradbuf[k], 0);
    }
    return 0;
};
fatrop_int OCPAdapter::eval_obj(
    OCPKKTMemory *OCP,
    double obj_scale,
    const FatropVecBF &primal_vars,
    double &res)
{
    // horizon length
    fatrop_int K = OCP->K;
    // offsets
    const fatrop_int *offs_ux = (const fatrop_int *)OCP->aux.ux_offs.data();
    fatrop_int *offs_stageparams_p = (fatrop_int *)offs_stageparams.data();
    double *stageparams_p = (double *)stageparams.data();
    double *globalparams_p = (double *)globalparams.data();
    OCPMACRO(fatrop_int *, nu, _p);
    SOLVERMACRO(VEC *, primal_vars, _p);
    double *primal_data = primal_vars_p->pa;
    double restot = 0.0;
    for (fatrop_int k = 0; k < K; k++)
    {
        fatrop_int nu_k = nu_p[k];
        fatrop_int offs_ux_k = offs_ux[k];
        fatrop_int offs_stageparams_k = offs_stageparams_p[k];
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

fatrop_int OCPAdapter::integrate_dynamics(
    OCPKKTMemory *OCP,
    const fatrop_int k,
    const FatropVecBF &uk,
    const FatropVecBF &xk,
    FatropVecBF &xkp1)
{
    // offsets
    fatrop_int *offs_stageparams_p = (fatrop_int *)offs_stageparams.data();
    fatrop_int offs_stageparams_k = offs_stageparams_p[k];
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
void OCPAdapter::set_parameters(const vector<double> &stage_params_in, const vector<double> &global_params_in)
{
    stageparams = stage_params_in;
    globalparams = global_params_in;
    return;
}
void OCPAdapter::set_initial_sol_guess(const shared_ptr<FatropData> &fatropdata, vector<double> &initial_u, vector<double> &initial_x)
{
    // offsets
    VEC *ux_intial_p = (VEC *)fatropdata->x_initial;
    double *u_p = initial_u.data();
    double *x_p = initial_x.data();
    fatrop_int offs_nu = 0;
    fatrop_int offs_nx = 0;
    fatrop_int offs_nux = 0;
    for (fatrop_int k = 0; k < K; k++)
    {
        fatrop_int nu_k = ocptempl->get_nuk(k);
        fatrop_int nx_k = ocptempl->get_nxk(k);
        PACKVEC(nu_k, u_p + offs_nu, 1, ux_intial_p, offs_nux);
        PACKVEC(nx_k, x_p + offs_nx, 1, ux_intial_p, offs_nux + nu_k);
        offs_nu += nu_k;
        offs_nx += nx_k;
        offs_nux += nu_k + nx_k;
    }
    return;
}
void OCPAdapter::get_solution(const shared_ptr<FatropData> &fatropdata, vector<double> &u, vector<double> &x)
{
    // offsets
    VEC *ux_sol = (VEC *)fatropdata->x_curr;
    double *u_p = u.data();
    double *x_p = x.data();
    fatrop_int offs_nu = 0;
    fatrop_int offs_nx = 0;
    fatrop_int offs_nux = 0;
    for (fatrop_int k = 0; k < K; k++)
    {
        fatrop_int nu_k = ocptempl->get_nuk(k);
        fatrop_int nx_k = ocptempl->get_nxk(k);
        UNPACKVEC(nu_k, ux_sol, offs_nux, u_p + offs_nu, 1);
        UNPACKVEC(nx_k, ux_sol, offs_nux + nu_k, x_p + offs_nx, 1);
        offs_nu += nu_k;
        offs_nx += nx_k;
        offs_nux += nu_k + nx_k;
    }
    return;
}