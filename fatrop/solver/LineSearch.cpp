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
    return fatropnlp_->EvalConstraintViolation(
        fatropdata_->x_next,
        fatropdata_->s_next,
        fatropdata_->g_next);
}
double LineSearch::EvalObjNext()
{
    double res = 0.0;
    fatropnlp_->EvalObj(
        fatropdata_->obj_scale,
        fatropdata_->x_next,
        res);
    return res;
}
int LineSearch::TryStep(double alpha_pr, double alpha_du) const
{
    return fatropdata_->TryStep(alpha_pr, alpha_du);
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
    fatropdata_->AlphaMax(alpha_primal, alpha_dual, MAX(1 - mu, 0.99));
    TryStep(alpha_primal, alpha_dual);
    double cv_curr = from_backup ? fatropdata_->CVL1Backup() : fatropdata_->CVL1Curr();
    double obj_curr = from_backup ? fatropdata_->obj_backup : fatropdata_->obj_curr;
    double barrier_curr = from_backup ? fatropdata_->EvalBarrierBackup(mu) : fatropdata_->EvalBarrierCurr(mu);
    obj_curr += barrier_curr;
    double lin_decr_curr = from_backup ? fatropdata_->LinDecrBackup() : fatropdata_->LinDecrCurr();
    double barrier_decr_curr = from_backup? fatropdata_->EvalBarrierLinDecrBackup(mu):  fatropdata_->EvalBarrierLinDecrCurr(mu);
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
        // todo change iteration number from zero to real iteration number
        if (small_sd)
        {
            (journaller_->it_curr).type = 's';
            fatropdata_->TakeStep();
            journaller_->it_curr.alpha_pr = alpha_primal;
            journaller_->it_curr.alpha_du = alpha_dual;
            res.ls = 1;
            return res;
        }
        if (filter_->IsAcceptable(FilterData(0, obj_next, cv_next)))
        {
            bool switch_cond = (lin_decr_curr < 0) && (alpha_primal * pow(-lin_decr_curr, s_phi) > delta * pow(cv_curr, s_theta));
            bool armijo = obj_next - obj_curr < eta_phi * alpha_primal * lin_decr_curr;
            if (switch_cond && (cv_curr <= theta_min))
            {
                // f-step
                if (armijo)
                {
                    (journaller_->it_curr).type = 'f';
                    fatropdata_->TakeStep();
                    journaller_->it_curr.alpha_pr = alpha_primal;
                    journaller_->it_curr.alpha_du = alpha_dual;
                    res.ls = ll;
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
                    (journaller_->it_curr).type = 'h';
                    fatropdata_->TakeStep();
                    journaller_->it_curr.alpha_pr = alpha_primal;
                    journaller_->it_curr.alpha_du = alpha_dual;
                    res.ls = ll;
                    return res;
                }
            }
        }
        else
        {
            if (ll == 1)
            {
                res.first_rejected_by_filter = true;
            }
        }
        alpha_primal *= 0.50;
    }
    assert(false);
    res.ls = 0;
    return res;
};