#include "BFOCPAL.hpp"
using namespace fatrop;
int BFOCPAL::get_nxk(const int k) const
{
    return ocp_->get_nxk(k);
};

int BFOCPAL::get_nuk(const int k) const
{
    return ocp_->get_nuk(k);
};

int BFOCPAL::get_ngk(const int k) const
{
    return ocp_->get_ngk(k);
}

int BFOCPAL::get_n_stage_params_k(const int k) const
{
    return ocp_->get_n_stage_params_k(k);
}

int BFOCPAL::get_n_global_parmas() const
{
    return ocp_->get_n_global_parmas();
}

int BFOCPAL::get_ng_ineq_k(const int k) const
{
    return 0;
}
int BFOCPAL::get_ng_ineq_k_AL(const int k) const
{
    return ocp_->get_ng_ineq_k(k);
}

int BFOCPAL::get_horizon_length() const
{
    return ocp_->get_horizon_length();
}

int BFOCPAL::eval_BAbtk(
    const double *states_kp1,
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params_k,
    MAT *res,
    const int k)
{
    ocp_->eval_BAbtk(
        states_kp1,
        inputs_k,
        states_k,
        stage_params_k,
        global_params_k,
        res,
        k);
    return 0;
};

int BFOCPAL::eval_RSQrqtk(
    const double *objective_scale,
    const double *inputs_k,
    const double *states_k,
    const double *lam_dyn_k,
    const double *lam_eq_k,
    const double *lam_eq_ineq_k,
    const double *stage_params_k,
    const double *global_params_k,
    MAT *res,
    const int k)
{
    ocp_->eval_RSQrqtk(
        objective_scale,
        inputs_k,
        states_k,
        lam_dyn_k,
        lam_eq_k,
        lam_eq_ineq_k,
        stage_params_k,
        global_params_k,
        res,
        k);
    // TODO add penalty
    return 0;
};

int BFOCPAL::eval_Ggtk(
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params_k,
    MAT *res,
    const int k)
{
    ocp_->eval_Ggtk(
        inputs_k,
        states_k,
        stage_params_k,
        global_params_k,
        res,
        k);
    return 0;
};

int BFOCPAL::eval_Ggt_ineqk(
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params_k,
    MAT *res,
    const int k)
{
    return 0;
};

int BFOCPAL::eval_bk(
    const double *states_kp1,
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params_k,
    double *res,
    const int k)
{
    ocp_->eval_bk(
        states_kp1,
        inputs_k,
        states_k,
        stage_params_k,
        global_params_k,
        res,
        k);
    return 0;
};

int BFOCPAL::eval_gk(
    const double *states_k,
    const double *inputs_k,
    const double *stage_params_k,
    const double *global_params_k,
    double *res,
    const int k)
{
    ocp_->eval_gk(
        states_k,
        inputs_k,
        stage_params_k,
        global_params_k,
        res,
        k);
    return 0;
};

int BFOCPAL::eval_gineqk(
    const double *states_k,
    const double *inputs_k,
    const double *stage_params_k,
    const double *global_params_k,
    double *res,
    const int k)
{
    return 0;
};
int BFOCPAL::eval_gineqk_AL(
    const double *states_k,
    const double *inputs_k,
    const double *stage_params_k,
    const double *global_params_k,
    double *res,
    const int k)
{
    ocp_->eval_gineqk(
        states_k,
        inputs_k,
        stage_params_k,
        global_params_k,
        res,
        k);
    return 0;
}

int BFOCPAL::eval_rqk(
    const double *objective_scale,
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params_k,
    double *res,
    const int k)
{
    ocp_->eval_rqk(
        objective_scale,
        inputs_k,
        states_k,
        stage_params_k,
        global_params_k,
        res,
        k);
    return 0;
};

int BFOCPAL::eval_Lk(
    const double *objective_scale,
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params_k,
    double *res,
    const int k)
{
    ocp_->eval_Lk(
        objective_scale,
        inputs_k,
        states_k,
        stage_params_k,
        global_params_k,
        res,
        k);
    // TODO add penalty
    return 0;
};