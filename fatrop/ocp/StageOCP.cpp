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
#include "ocp/StageOCP.hpp"
using namespace fatrop;
using namespace std;
StageOCPRockit::StageOCPRockit(
    const fatrop_int nu,
    const fatrop_int nx,
    const fatrop_int ngI,
    const fatrop_int ng,
    const fatrop_int ngF,
    const fatrop_int ng_ineqI,
    const fatrop_int ng_ineq,
    const fatrop_int ng_ineqF,
    const fatrop_int n_stage_params,
    const fatrop_int n_global_params,
    const fatrop_int K,
    const EvalCasGen &BAbtf,
    const EvalCasGen &bkf,
    const EvalCasGen &RSQrqtIf,
    const EvalCasGen &rqIf,
    const EvalCasGen &RSQrqtf,
    const EvalCasGen &rqf,
    const EvalCasGen &RSQrqtFf,
    const EvalCasGen &rqFf,
    const EvalCasGen &GgtIf,
    const EvalCasGen &gIf,
    const EvalCasGen &Ggtf,
    const EvalCasGen &gf,
    const EvalCasGen &GgtFf,
    const EvalCasGen &gFf,
    const EvalCasGen &Ggt_ineqIf,
    const EvalCasGen &gineqIf,
    const EvalCasGen &Ggt_ineqf,
    const EvalCasGen &gineqf,
    const EvalCasGen &Ggt_ineqFf,
    const EvalCasGen &gineqFf,
    const EvalCasGen &LkIf,
    const EvalCasGen &Lkf,
    const EvalCasGen &LFf,
    const vector<double> &bounds_L,
    const vector<double> &bounds_U,
    const vector<double> &stage_params,
    const vector<double> &global_params,
    const vector<double> &initial_u,
    const vector<double> &initial_x) : StageOCP(nu, nx, ngI, ng, ngF, ng_ineqI, ng_ineq, ng_ineqF, n_stage_params, n_global_params, K),
                                       BAbtf(BAbtf),
                                       bkf(bkf),
                                       RSQrqtIf(RSQrqtIf),
                                       rqIf(rqIf),
                                       RSQrqtf(RSQrqtf),
                                       rqf(rqf),
                                       RSQrqtFf(RSQrqtFf),
                                       rqFf(rqFf),
                                       GgtIf(GgtIf),
                                       gIf(gIf),
                                       Ggtf(Ggtf),
                                       gf(gf),
                                       GgtFf(GgtFf),
                                       gFf(gFf),
                                       Ggt_ineqIf(Ggt_ineqIf),
                                       g_ineqIf(gineqIf),
                                       Ggt_ineqf(Ggt_ineqf),
                                       g_ineqf(gineqf),
                                       Ggt_ineqFf(Ggt_ineqFf),
                                       g_ineqFf(gineqFf),
                                       LkIf(LkIf),
                                       Lkf(Lkf),
                                       LFf(LFf),
                                       initial_x(initial_x),
                                       initial_u(initial_u),
                                       bounds_L(bounds_L),
                                       bounds_U(bounds_U),
                                       stage_params(stage_params),
                                       global_params(global_params)
{
}
fatrop_int StageOCP::get_nxk(const fatrop_int k) const
{
    return nx_;
}
fatrop_int StageOCP::get_nuk(const fatrop_int k) const
{
    if (k == K_ - 1)
        return 0;
    return nu_;
}
fatrop_int StageOCP::get_ngk(const fatrop_int k) const
{
    if (k == 0)
        return ngI_;
    if (k == K_ - 1)
        return ngF_;
    return ng_;
}
fatrop_int StageOCP::get_ng_ineq_k(const fatrop_int k) const
{
    if (k == 0)
    {
        return ng_ineqI_;
    }
    if (k == K_ - 1)
    {
        return ng_ineqF_;
    }
    return ng_ineq_;
}
fatrop_int StageOCP::get_n_global_params() const
{
    return n_global_params_;
};
fatrop_int StageOCP::get_n_stage_params_k(const fatrop_int k) const
{
    return n_stage_params_;
};
fatrop_int StageOCP::get_horizon_length() const
{
    return K_;
};
fatrop_int StageOCPRockit::eval_BAbtk(const double *states_kp1,
                               const double *inputs_k,
                               const double *states_k,
                               const double *stage_params_k,
                               const double *global_params,
                               MAT *res,
                               const fatrop_int k)
{
    const double *args[5];
    args[0] = states_kp1;
    args[1] = inputs_k;
    args[2] = states_k;
    args[3] = stage_params_k;
    args[4] = global_params;
    return BAbtf.eval_bf(args, res);
}
fatrop_int StageOCPRockit::eval_RSQrqtk(const double *objective_scale,
                                 const double *inputs_k,
                                 const double *states_k,
                                 const double *lam_dyn_k,
                                 const double *lam_eq_k,
                                 const double *lam_ineq_k,
                                 const double *stage_params_k,
                                 const double *global_params,
                                 MAT *res,
                                 const fatrop_int k)
{
    const double *args[8];
    args[0] = objective_scale;
    args[1] = inputs_k;
    args[2] = states_k;
    args[3] = lam_dyn_k;
    args[4] = lam_eq_k;
    args[5] = lam_ineq_k;
    args[6] = stage_params_k;
    args[7] = global_params;
    if (k == 0)
        return RSQrqtIf.eval_bf(args, res);
    if (k == K_ - 1)
        return RSQrqtFf.eval_bf(args, res);
    return RSQrqtf.eval_bf(args, res);
};
fatrop_int StageOCPRockit::eval_Ggtk(
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params,
    MAT *res,
    const fatrop_int k)
{
    const double *args[4];
    args[0] = inputs_k;
    args[1] = states_k;
    args[2] = stage_params_k;
    args[3] = global_params;
    if (k == K_ - 1)
        return GgtFf.eval_bf(args, res);
    if (k == 0)
        return GgtIf.eval_bf(args, res);
    return Ggtf.eval_bf(args, res);
};
fatrop_int StageOCPRockit::eval_Ggt_ineqk(
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params,
    MAT *res,
    const fatrop_int k)
{
    if (k == K_ - 1)
    {
        const double *args[3];
        args[0] = states_k;
        args[1] = stage_params_k;
        args[2] = global_params;
        return Ggt_ineqFf.eval_bf(args, res);
    }
    const double *args[4];
    args[0] = inputs_k;
    args[1] = states_k;
    args[2] = stage_params_k;
    args[3] = global_params;
    if (k == 0)
        return Ggt_ineqIf.eval_bf(args, res);
    return Ggt_ineqf.eval_bf(args, res);
};
fatrop_int StageOCPRockit::eval_bk(
    const double *states_kp1,
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params,
    double *constraint_violation_k,
    const fatrop_int k)
{
    const double *args[5];
    args[0] = states_kp1;
    args[1] = inputs_k;
    args[2] = states_k;
    args[3] = stage_params_k;
    args[4] = global_params;
    return bkf.eval_array(args, constraint_violation_k);
};
fatrop_int StageOCPRockit::eval_gk(
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params,
    double *res,
    const fatrop_int k)
{
    const double *args[4];
    args[0] = inputs_k;
    args[1] = states_k;
    args[2] = stage_params_k;
    args[3] = global_params;
    if (k == K_ - 1)
        return gFf.eval_array(args, res);
    if (k == 0)
        return gIf.eval_array(args, res);
    return gf.eval_array(args, res);
}
fatrop_int StageOCPRockit::eval_gineqk(
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params,
    double *res,
    const fatrop_int k)
{
    if (k == K_ - 1)
    {
        const double *args[3];
        args[0] = states_k;
        args[1] = stage_params_k;
        args[2] = global_params;
        return g_ineqFf.eval_array(args, res);
    }
    const double *args[4];
    args[0] = inputs_k;
    args[1] = states_k;
    args[2] = stage_params_k;
    args[3] = global_params;
    if (k == 0)
        return g_ineqIf.eval_array(args, res);
    return g_ineqf.eval_array(args, res);
}
fatrop_int StageOCPRockit::eval_rqk(
    const double *objective_scale,
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params,
    double *res,
    const fatrop_int k)
{
    const double *args[5];
    args[0] = objective_scale;
    args[1] = inputs_k;
    args[2] = states_k;
    args[3] = stage_params_k;
    args[4] = global_params;
    if (k == K_ - 1)
        return rqFf.eval_array(args, res);
    if (k == 0)
        return rqIf.eval_array(args, res);
    return rqf.eval_array(args, res);
};
fatrop_int StageOCPRockit::eval_Lk(
    const double *objective_scale,
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params,
    double *res,
    const fatrop_int k)
{
    const double *args[5];
    args[0] = objective_scale;
    args[1] = inputs_k;
    args[2] = states_k;
    args[3] = stage_params_k;
    args[4] = global_params;
    if (k == K_ - 1)
        return LFf.eval_array(args, res);
    if (k == 0)
        return LkIf.eval_array(args, res);
    return Lkf.eval_array(args, res);
};