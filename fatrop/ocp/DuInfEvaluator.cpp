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
#include "ocp/DuInfEvaluator.hpp"
using namespace fatrop;
fatrop_int DuInfEvaluator::evaluate(
    OCPKKTMemory *OCP,
    double obj_scale,
    const FatropVecBF &lam,
    const FatropVecBF &grad,
    FatropVecBF &du_inf)
{
    fatrop_int K = OCP->K;
    VEC *lam_p = (VEC *)lam;
    VEC *grad_p = (VEC *)grad;
    VEC *du_inf_p = (VEC *)du_inf;
    OCPMACRO(MAT *, BAbt, _p);
    OCPMACRO(MAT *, Ggt, _p);
    OCPMACRO(MAT *, Ggt_ineq, _p);
    OCPMACRO(fatrop_int *, nu, _p);
    OCPMACRO(fatrop_int *, nx, _p);
    OCPMACRO(fatrop_int *, ng, _p);
    OCPMACRO(fatrop_int *, ng_ineq, _p);
    DBGASSERT(grad_p->m == du_inf_p->m);
    // copy grad_f to du_inf
    VECCP(grad_p->m, grad_p, 0, du_inf_p, 0);
    fatrop_int *offs_ux = (fatrop_int *)OCP->aux.ux_offs.data();
    fatrop_int *offs_g = (fatrop_int *)OCP->aux.g_offs.data();
    fatrop_int *offs_dyn_eq = (fatrop_int *)OCP->aux.dyn_eq_offs.data();
    fatrop_int *offs_g_ineq = (fatrop_int *)OCP->aux.g_ineq_offs.data();
    // contribution of dynamics constraints
    for (fatrop_int k = 0; k < K - 1; k++)
    {
        const fatrop_int nu = nu_p[k];
        const fatrop_int nup1 = nu_p[k + 1];
        const fatrop_int nx = nx_p[k];
        const fatrop_int nxp1 = nx_p[k + 1];
        const fatrop_int offsp1 = offs_ux[k + 1];
        const fatrop_int offs = offs_ux[k];
        const fatrop_int offs_dyn_eq_k = offs_dyn_eq[k];
        GEMV_N(nu + nx, nxp1, 1.0, BAbt_p + k, 0, 0, lam_p, offs_dyn_eq_k, 1.0, du_inf_p, offs, du_inf_p, offs);
        AXPY(nxp1, -1.0, lam_p, offs_dyn_eq_k, du_inf_p, offsp1 + nup1, du_inf_p, offsp1 + nup1);
    }
    // contribution of equality constraints
    for (fatrop_int k = 0; k < K; k++)
    {
        const fatrop_int nu = nu_p[k];
        const fatrop_int nx = nx_p[k];
        const fatrop_int ng = ng_p[k];
        const fatrop_int offs = offs_ux[k];
        const fatrop_int offs_g_k = offs_g[k];
        GEMV_N(nu + nx, ng, 1.0, Ggt_p + k, 0, 0, lam_p, offs_g_k, 1.0, du_inf_p, offs, du_inf_p, offs);
    }
    // constribution of inequality - slack constraints
    for (fatrop_int k = 0; k < K; k++)
    {
        const fatrop_int nu = nu_p[k];
        const fatrop_int nx = nx_p[k];
        const fatrop_int ng_ineq = ng_ineq_p[k];
        const fatrop_int offs_g_ineq_k = offs_g_ineq[k];
        const fatrop_int offs = offs_ux[k];
        GEMV_N(nu + nx, ng_ineq, 1.0, Ggt_ineq_p + k, 0, 0, lam_p, offs_g_ineq_k, 1.0, du_inf_p, offs, du_inf_p, offs);
    }
    return 0;
}