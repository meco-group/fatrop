#include "solver/LineSearch.hpp"
using namespace fatrop;
LineSearch::LineSearch(
    const shared_ptr<FatropParams> &fatropparams,
    const shared_ptr<FatropNLP> &nlp,
    const shared_ptr<FatropData> &fatropdata) : AlgStrategy(fatropparams),
                                                fatropnlp_(nlp),
                                                fatropdata_(fatropdata){};
inline int LineSearch::EvalCVNext()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    int res = fatropnlp_->EvalConstraintViolation(
        fatropdata_->x_next,
        fatropdata_->s_next,
        fatropdata_->g_next);
    eval_cv_time += blasfeo_toc(&timer);
    eval_cv_count++;
    return res;
}
double LineSearch::EvalObjNext()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    double res = 0.0;
    fatropnlp_->EvalObj(
        fatropdata_->obj_scale,
        fatropdata_->x_next,
        res);
    eval_obj_time += blasfeo_toc(&timer);
    eval_obj_count++;
    return res;
}
void LineSearch::Reset()
{
    eval_cv_count = 0;
    eval_obj_count = 0;
    eval_cv_time = 0.;
    eval_obj_time = 0.;
}
int LineSearch::TryStep(double alpha_pr, double alpha_du) const
{
    return fatropdata_->TryStep(alpha_pr, alpha_du);
};
int LineSearch::InitSoc() const
{
    // backup delta_x, delta_s and lam_calc
    fatropdata_->lam_calc_backup_ls.copy(fatropdata_->lam_calc);
    fatropdata_->delta_x_backup_ls.copy(fatropdata_->delta_x);
    fatropdata_->delta_s_backup_ls.copy(fatropdata_->delta_s);
    fatropdata_->g_soc.copy(fatropdata_->g_curr);
    return 0;
};
int LineSearch::ExitSoc() const
{
    // restore delta_x, delta_s and lam_calc
    fatropdata_->lam_calc.copy(fatropdata_->lam_calc_backup_ls);
    fatropdata_->delta_x.copy(fatropdata_->delta_x_backup_ls);
    fatropdata_->delta_s.copy(fatropdata_->delta_s_backup_ls);
    return 0;
};
int LineSearch::CalcSoc(double alpha) const
{
    axpy(alpha, fatropdata_->g_soc, fatropdata_->g_next, fatropdata_->g_soc);
    return fatropnlp_->SolveSOC(fatropdata_->delta_x, fatropdata_->lam_calc, fatropdata_->delta_s, fatropdata_->g_soc);
};

BackTrackingLineSearch::BackTrackingLineSearch(
    const shared_ptr<FatropParams> &fatropparams,
    const shared_ptr<FatropNLP> &nlp,
    const shared_ptr<FatropData> &fatropdata,
    const shared_ptr<Filter> &filter,
    const shared_ptr<Journaller> &journaller)
    : LineSearch(fatropparams, nlp, fatropdata), filter_(filter), journaller_(journaller)
{
    Initialize();
};
void BackTrackingLineSearch::Initialize()
{
    // todo avoid reallocation when maxiter doesn't change
    s_phi = fatrop_params_->s_phi;
    delta = fatrop_params_->delta;
    s_theta = fatrop_params_->s_theta;
    gamma_theta = fatrop_params_->gamma_theta;
    gamma_phi = fatrop_params_->gamma_phi;
    eta_phi = fatrop_params_->eta_phi;
    gamma_alpha = fatrop_params_->gamma_alpha;
}
LineSearchInfo BackTrackingLineSearch::FindAcceptableTrialPoint(double mu, bool small_sd, bool from_backup)
{
    LineSearchInfo res;
    double alpha_primal = 1.0;
    double alpha_dual = 1.0;
    double alpha_primal_backup = 1.0;
    double alpha_dual_backup = 1.0;
    fatropdata_->AlphaMax(alpha_primal, alpha_dual, MAX(1 - mu, 0.99));
    TryStep(alpha_primal, alpha_dual);
    double cv_curr = from_backup ? fatropdata_->CVL1Backup() : fatropdata_->CVL1Curr();
    double obj_curr = from_backup ? fatropdata_->obj_backup : fatropdata_->obj_curr;
    double barrier_curr = from_backup ? fatropdata_->EvalBarrierBackup(mu) : fatropdata_->EvalBarrierCurr(mu);
    obj_curr += barrier_curr;
    double lin_decr_curr = from_backup ? fatropdata_->LinDecrBackup() : fatropdata_->LinDecrCurr();
    double barrier_decr_curr = from_backup ? fatropdata_->EvalBarrierLinDecrBackup(mu) : fatropdata_->EvalBarrierLinDecrCurr(mu);
    lin_decr_curr += barrier_decr_curr;
    double theta_min = fatropdata_->theta_min;
    // calculation of alpha_min
    double alpha_min = gamma_alpha *
                       (lin_decr_curr > 0 ? gamma_theta
                                          : MIN(gamma_theta,
                                                cv_curr < theta_min ? MIN(-gamma_phi * cv_curr / lin_decr_curr, delta * pow(cv_curr, s_theta) / (pow(-lin_decr_curr, s_phi)))
                                                                    : -gamma_phi * cv_curr / lin_decr_curr));

    // cout << "alpha_min " << alpha_min << endl;
    // cout << "cv " << cv_curr << endl;
    // cout << "obj " << obj_curr << endl;
    // cout << "lindecr " << lin_decr_curr << endl;
    const int p_max = 0;
    bool soc_step = false;
    double cv_soc_old = cv_curr;
    int p = 0;
    int no_alpha_trials = 1;
    for (int ll = 1; ll < 50; ll++)
    {
        TryStep(alpha_primal, alpha_dual);
        if (alpha_primal < alpha_min)
        {
            res.ls = 0;
            return res;
        }
        EvalCVNext();
        double cv_next = fatropdata_->CVL1Next();
        double obj_next = EvalObjNext();
        double barrier_next = fatropdata_->EvalBarrierNext(mu);
        obj_next += barrier_next;

// cout << "cv_next: " << cv_next;
// cout << "  obj_next: " << obj_next << endl;
        double alpha_primal_accent = (soc_step ? alpha_primal_backup : alpha_primal);
        if (filter_->IsAcceptable(FilterData(0, obj_next, cv_next)))
        {
            // cout << filter_->GetSize() << endl;
            bool switch_cond = (lin_decr_curr < 0) && (alpha_primal_accent* pow(-lin_decr_curr, s_phi) > delta * pow(cv_curr, s_theta));
            bool armijo = obj_next - obj_curr < eta_phi *alpha_primal_accent * lin_decr_curr;
            if (switch_cond && (cv_curr <= theta_min))
            {
                // f-step
                if (armijo)
                {
                    (journaller_->it_curr).type = soc_step? 'F':'f';
                    fatropdata_->TakeStep();
                    journaller_->it_curr.alpha_pr = alpha_primal;
                    journaller_->it_curr.alpha_du = alpha_dual;
                    res.ls = no_alpha_trials;
                    return res;
                }
            }
            else
            {
                // h-step
                // check sufficient decrease wrt current iterate
                if ((cv_next < (1.0 - gamma_theta) * cv_curr) || (obj_next < obj_curr - gamma_phi * cv_curr))
                {
                    if (!switch_cond || !(armijo))
                    {
                        filter_->Augment(FilterData(0, obj_curr - gamma_phi * cv_curr, cv_curr * (1 - gamma_theta)));
                    }
                    (journaller_->it_curr).type = soc_step? 'H':'h';
                    fatropdata_->TakeStep();
                    journaller_->it_curr.alpha_pr = alpha_primal;
                    journaller_->it_curr.alpha_du = alpha_dual;
                    res.ls = no_alpha_trials;
                    return res;
                }
            }
        }
        else
        {
            if(soc_step)
            {
                // abort soc
                p = p_max;
            }
            if (ll == 1)
            {
                res.first_rejected_by_filter = true;
            }
        }
        // todo change iteration number from zero to real iteration number
        if (small_sd)
        {
            (journaller_->it_curr).type = 's';
            fatropdata_->TakeStep();
            journaller_->it_curr.alpha_pr = alpha_primal;
            journaller_->it_curr.alpha_du = alpha_dual;
            res.ls = -1;
            return res;
        }
        if(soc_step && (p>p_max || (ll>1 &&(cv_next > 0.99*cv_soc_old))))
        {
            // deactivate soc
            soc_step = false;
            alpha_primal = alpha_primal_backup;
            alpha_dual = alpha_dual_backup;
            ExitSoc();
            // todo cache these variables
            fatropdata_->ComputedZ();
        }
        if (!soc_step && (ll == 1 && p_max>0))
        {
            // activate soc
            cout << "trying soc " << endl;
            soc_step = true;
            alpha_primal_backup = alpha_primal;
            alpha_dual_backup = alpha_dual;
            InitSoc();
        }
        if(soc_step)
        {
            if(ll>1)
            {
                cv_soc_old = cv_next;
            }
            CalcSoc(alpha_primal);
            // cout << "size of soc x step " << L1(fatropdata_->delta_x) << endl;
            fatropdata_->ComputedZ();
            fatropdata_->AlphaMax(alpha_primal, alpha_dual, MAX(1 - mu, 0.99));
            // cout << "alpha_primal " << alpha_primal << endl;
            // cout << "alpha_dual " << alpha_dual << endl;
            // cout << "soc" << endl;
            p = p + 1;
        }
        else
        {
            // cut back alpha_primal
            no_alpha_trials++;
            alpha_primal *= 0.50;
        }

    }
    assert(false);
    res.ls = 0;
    return res;
};