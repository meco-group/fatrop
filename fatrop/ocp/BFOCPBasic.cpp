#include "ocp/BFOCPBasic.hpp"
using namespace fatrop;
BFOCPBasic::BFOCPBasic(const int nu,
                       const int nx,
                       const int ngI,
                       const int ng,
                       const int ngF,
                       const int ng_ineqI,
                       const int ng_ineq,
                       const int ng_ineqF,
                       const int n_stage_params,
                       const int n_global_params,
                       const int K,
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
                       const EvalCasGen &LFf) : nu_(nu),
                                                nx_(nx),
                                                ngI_(ngI),
                                                ng_(ng),
                                                ngF_(ngF),
                                                ng_ineqI_(ng_ineqI),
                                                ng_ineq_(ng_ineq),
                                                ng_ineqF_(ng_ineqF),
                                                n_stage_params_(n_stage_params),
                                                n_global_params_(n_global_params),
                                                K_(K),
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
                                                initial_x(K * nx, 0.0),
                                                initial_u((K - 1) * nu_, 0.0)
{
}
int BFOCPBasic::get_nxk(const int k) const
{
    return nx_;
}
int BFOCPBasic::get_nuk(const int k) const
{
    if (k == K_ - 1)
        return 0;
    return nu_;
}
int BFOCPBasic::get_ngk(const int k) const
{
    if (k == 0)
        return ngI_;
    if (k == K_ - 1)
        return ngF_;
    return ng_;
}
int BFOCPBasic::get_ng_ineq_k(const int k) const
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
int BFOCPBasic::get_n_global_params() const
{
    return n_global_params_;
};
int BFOCPBasic::get_n_stage_params_k(const int k) const
{
    return n_stage_params_;
};
int BFOCPBasic::get_horizon_length() const
{
    return K_;
};
int BFOCPBasic::eval_BAbtk(const double *states_kp1,
                           const double *inputs_k,
                           const double *states_k,
                           const double *stage_params_k,
                           const double *global_params,
                           MAT *res,
                           const int k)
{
    const double *args[5];
    args[0] = states_kp1;
    args[1] = inputs_k;
    args[2] = states_k;
    args[3] = stage_params_k;
    args[4] = global_params;
    return BAbtf.eval_bf(args, res);
}
int BFOCPBasic::eval_RSQrqtk(const double *objective_scale,
                             const double *inputs_k,
                             const double *states_k,
                             const double *lam_dyn_k,
                             const double *lam_eq_k,
                             const double *lam_ineq_k,
                             const double *stage_params_k,
                             const double *global_params,
                             MAT *res,
                             const int k)
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
int BFOCPBasic::eval_Ggtk(
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params,
    MAT *res,
    const int k)
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
int BFOCPBasic::eval_Ggt_ineqk(
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params,
    MAT *res,
    const int k)
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
int BFOCPBasic::eval_bk(
    const double *states_kp1,
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params,
    double *constraint_violation_k,
    const int k)
{
    const double *args[5];
    args[0] = states_kp1;
    args[1] = inputs_k;
    args[2] = states_k;
    args[3] = stage_params_k;
    args[4] = global_params;
    return bkf.eval_array(args, constraint_violation_k);
};
int BFOCPBasic::eval_gk(
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params,
    double *res,
    const int k)
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
int BFOCPBasic::eval_gineqk(
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params,
    double *res,
    const int k)
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
int BFOCPBasic::eval_rqk(
    const double *objective_scale,
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params,
    double *res,
    const int k)
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
int BFOCPBasic::eval_Lk(
    const double *objective_scale,
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params,
    double *res,
    const int k)
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