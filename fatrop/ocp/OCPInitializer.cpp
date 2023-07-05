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
#include "ocp/OCPInitializer.hpp"
using namespace fatrop;

fatrop_int OCPInitializer::modify_kkt_ls_dual_estimate(
    OCPKKTMemory *OCP,
    const FatropVecBF &grad)
{
    // horizon length
    fatrop_int K = OCP->K;
    // offsets
    fatrop_int *offs_ux = (fatrop_int *)OCP->aux.ux_offs.data();
    OCPMACRO(MAT *, BAbt, _p);
    OCPMACRO(MAT *, Ggt, _p);
    OCPMACRO(MAT *, RSQrqt, _p);
    OCPMACRO(fatrop_int *, nu, _p);
    OCPMACRO(fatrop_int *, nx, _p);
    OCPMACRO(fatrop_int *, ng, _p);
    SOLVERMACRO(VEC *, grad, _p);
    // SOLVERMACRO(VEC *, s, _p);
    for (fatrop_int k = 0; k < K; k++)
    {
        fatrop_int nu_k = nu_p[k];
        fatrop_int nx_k = nx_p[k];
        fatrop_int offs_ux_k = offs_ux[k];
        fatrop_identity(nu_k + nx_k, RSQrqt_p + k, 0, 0);
        ROWIN(nu_k + nx_k, 1.0, grad_p, offs_ux_k, RSQrqt_p + k, nu_k + nx_k, 0);
    }

    for (fatrop_int k = 0; k < K - 1; k++)
    {
        fatrop_int nu_k = nu_p[k];
        fatrop_int nx_k = nx_p[k];
        fatrop_int nx_kp1 = nx_p[k + 1];
        GESE(1, nx_kp1, 0.0, BAbt_p + k, nu_k + nx_k, 0);
    }
    for (fatrop_int k = 0; k < K; k++)
    {
        fatrop_int nu_k = nu_p[k];
        fatrop_int nx_k = nx_p[k];
        fatrop_int ng_k = ng_p[k];
        if (ng_k > 0)
        {
            GESE(1, ng_k, 0.0, Ggt_p + k, nu_k + nx_k, 0);
        }
    }
    OCPMACRO(MAT *, Ggt_ineq, _p);
    OCPMACRO(fatrop_int *, ng_ineq, _p);
    for (fatrop_int k = 0; k < K; k++)
    {
        fatrop_int nu_k = nu_p[k];
        fatrop_int nx_k = nx_p[k];
        fatrop_int ng_ineq_k = ng_ineq_p[k];
        if (ng_ineq_k > 0)
        {
            GESE(1, ng_ineq_k, 0.0, Ggt_ineq_p + k, nu_k + nx_k, 0);
        }
    }
    return 0;
}
fatrop_int OCPInitializer::intialize_slack_variables(
    OCPKKTMemory *OCP,
    FatropVecBF &s)
{
    // horizon length
    fatrop_int K = OCP->K;
    // offsets
    fatrop_int *offs_ineq_p = (fatrop_int *)OCP->aux.ineq_offs.data();
    OCPMACRO(MAT *, Ggt_ineq, _p);
    OCPMACRO(fatrop_int *, nu, _p);
    OCPMACRO(fatrop_int *, nx, _p);
    OCPMACRO(fatrop_int *, ng_ineq, _p);
    SOLVERMACRO(VEC *, s, _p);
    for (fatrop_int k = 0; k < K; k++)
    {
        fatrop_int nu_k = nu_p[k];
        fatrop_int nx_k = nx_p[k];
        fatrop_int ng_ineq_k = ng_ineq_p[k];
        fatrop_int offs_ineq_k = offs_ineq_p[k];
        if (ng_ineq_k > 0)
        {
            ROWEX(ng_ineq_k, 1.0, Ggt_ineq_p + k, nu_k + nx_k, 0, s_p, offs_ineq_k);
            // GESE(1, ng_ineq_k, 0.0, Ggt_ineq_p + k, nu_k + nx_k, 0);
        }
    }
    return 0;
}