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
    int no_ineqsk = no_ineqs.at(k);
    int nuk = nu.at(k);
    int nxk = nx.at(k);
    int offs = ineqs_offsets.at(k);
    double *lowerp =lower_bounds.data();
    double *upperp =upper_bounds.data();
    double *tmpviolationp = ((blasfeo_dvec *)tmpviolation)->pa;
    double *ineq_lagsLp = ((blasfeo_dvec *)ineq_lagsL)->pa + offs;
    double *ineq_lagsUp = ((blasfeo_dvec *)ineq_lagsU)->pa + offs;
    double *lagsupdatedp = ((blasfeo_dvec *)lagsupdated)->pa;
    double penalty = this->penalty;
    MAT *tmpmatp = (MAT *)this->tmpmat;
    // evaluate ineq Jacobian
    this->eval_Ggt_ineqk_AL(
        inputs_k,
        states_k,
        stage_params_k,
        global_params_k,
        tmpmatp,
        k);
    // evaluate ineq vector
    // todo with ROWEX
    this->eval_gineqk_AL(
        states_k,
        inputs_k,
        stage_params_k,
        global_params_k,
        tmpviolationp,
        k);
    // calculate updated lagineqs and vector for gradient
    for (int i = 0; i < no_ineqsk; i++)
    {
        double ineqlagsLi = ineq_lagsLp[i];
        double ineqlagsUi = ineq_lagsUp[i];
        // lambdai = 0
        double violationi = tmpviolationp[i];
        double loweri = lowerp[i];
        bool lower_bounded = !isinf(loweri);
        double upperi = upperp[i];
        bool upper_bounded = !isinf(upperi);
        double dist_low = lower_bounded ? -ineqlagsLi + penalty*(violationi - loweri) : 0.0;
        double dist_up = upper_bounded ? -ineqlagsUi +  penalty*(upperi - violationi): 0.0;
        lagsupdatedp[i] = 0;
        if (dist_low < 0.0)
        {
            lagsupdatedp[i] += dist_low;
        }
        if (dist_up < 0.0)
        {
            lagsupdatedp[i] += -dist_up;
        }
        if (dist_low < 0.0 && dist_up < 0.0)
        {
            GESC(nuk + nxk, 1, sqrt(2), tmpmatp, 0, 0);
        }
        if (dist_low >= 0.0 && dist_up >= 0.0)
        {
            // turn off penalty
            GESE(nuk + nxk, 1, 0.0, tmpmatp, 0, 0);
        }
        // gradvecp[i] =;
    }
    ocp_->eval_RSQrqtk(
        objective_scale,
        inputs_k,
        states_k,
        lam_dyn_k,
        lam_eq_k,
        lagsupdatedp,
        stage_params_k,
        global_params_k,
        res,
        k);
    SYRK_LN_MN(nuk + nxk, nuk + nxk, no_ineqsk, penalty, tmpmatp, 0, 0, tmpmatp, 0, 0, 1.0, res, 0, 0, res, 0, 0);
    TRTR_L(nuk + nxk, res, 0, 0, res, 0, 0);
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
int BFOCPAL::eval_Ggt_ineqk_AL(
    const double *inputs_k,
    const double *states_k,
    const double *stage_params_k,
    const double *global_params_k,
    MAT *res,
    const int k)
{
    ocp_->eval_Ggt_ineqk(
        inputs_k,
        states_k,
        stage_params_k,
        global_params_k,
        res,
        k);
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
    int no_ineqsk = no_ineqs.at(k);
    int nuk = nu.at(k);
    int nxk = nx.at(k);
    int offs = ineqs_offsets.at(k);
    double *lowerp = lower_bounds.data();
    double *upperp = upper_bounds.data();
    double *tmpviolationp = ((blasfeo_dvec *)tmpviolation)->pa;
    double *ineq_lagsLp = ((blasfeo_dvec *)ineq_lagsL)->pa + offs;
    double *ineq_lagsUp = ((blasfeo_dvec *)ineq_lagsU)->pa + offs;
    double *gradvecp = ((blasfeo_dvec *)gradvec)->pa;
    double penalty = this->penalty;
    MAT *tmpmatp = (MAT *)this->tmpmat;
    // evaluate ineq Jacobian
    this->eval_Ggt_ineqk_AL(
        inputs_k,
        states_k,
        stage_params_k,
        global_params_k,
        tmpmatp,
        k);
    // evaluate ineq vector
    // todo with ROWEX
    this->eval_gineqk_AL(
        states_k,
        inputs_k,
        stage_params_k,
        global_params_k,
        tmpviolationp,
        k);
    // calculate updated lagineqs and vector for gradient
    for (int i = 0; i < no_ineqsk; i++)
    {
        double ineqlagsLi = ineq_lagsLp[i];
        double ineqlagsUi = ineq_lagsUp[i];
        // lambdai = 0
        double violationi = tmpviolationp[i];
        double loweri = lowerp[i];
        bool lower_bounded = !isinf(loweri);
        double upperi = upperp[i];
        bool upper_bounded = !isinf(upperi);
        double dist_low = lower_bounded ? -ineqlagsLi + penalty*(violationi - loweri) : 0.0;
        double dist_up = upper_bounded ? -ineqlagsUi +  penalty*(upperi - violationi): 0.0;
        gradvecp[i] = 0;
        if (dist_low < 0.0)
        {
            gradvecp[i] += dist_low;
        }
        if (dist_up < 0.0)
        {
            gradvecp[i] += -dist_up;
        }
    }
    ocp_->eval_rqk(
        objective_scale,
        inputs_k,
        states_k,
        stage_params_k,
        global_params_k,
        res,
        k);
    VEC resbf;
    CREATE_VEC(no_ineqsk, &resbf, res);
    // TODO add penalty
    GEMV_N(nuk + nxk, no_ineqsk, 1.0, tmpmatp, 0, 0, (VEC *)gradvec, 0, 1.0, &resbf, 0, &resbf, 0);
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
    int no_ineqsk = no_ineqs.at(k);
    int offs = ineqs_offsets.at(k);
    double *lowerp = lower_bounds.data();
    double *upperp = upper_bounds.data();
    double *tmpviolationp = ((blasfeo_dvec *)tmpviolation)->pa;
    double *ineq_lagsLp = ((blasfeo_dvec *)ineq_lagsL)->pa + offs;
    double *ineq_lagsUp = ((blasfeo_dvec *)ineq_lagsU)->pa + offs;
    double penalty = this->penalty;
    double penaltym1 = 1.0/penalty;
    double obj_penalty = 0.0;
    this->eval_gineqk_AL(
        states_k,
        inputs_k,
        stage_params_k,
        global_params_k,
        tmpviolationp,
        k);
    // calculate updated lagineqs and vector for gradient
    for (int i = 0; i < no_ineqsk; i++)
    {
        double ineqlagsLi = ineq_lagsLp[i];
        double ineqlagsUi = ineq_lagsUp[i];
        // lambdai = 0
        double violationi = tmpviolationp[i];
        double loweri = lowerp[i];
        bool lower_bounded = !isinf(loweri);
        double upperi = upperp[i];
        bool upper_bounded = !isinf(upperi);
        double dist_low = lower_bounded ? -ineqlagsLi + penalty*(violationi - loweri) : 0.0;
        double dist_up = upper_bounded ? -ineqlagsUi +  penalty*(upperi - violationi): 0.0;
        if (dist_low < 0.0)
        {
            obj_penalty += 0.5*penaltym1*dist_low*dist_low;
        }
        if (dist_up < 0.0)
        {
            obj_penalty += 0.5*penaltym1*dist_up*dist_up;
        }
    }
    ocp_->eval_Lk(
        objective_scale,
        inputs_k,
        states_k,
        stage_params_k,
        global_params_k,
        res,
        k);
    *res += obj_penalty;
    return 0;
};