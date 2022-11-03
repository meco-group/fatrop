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
int BackTrackingLineSearch::FindAcceptableTrialPoint(double mu, bool small_sd)
{
    double alpha_primal = 1.0;
    double alpha_dual = 1.0;
    fatropdata_->AlphaMax(alpha_primal, alpha_dual, MAX(1 - mu, 0.99));
    TryStep(alpha_primal, alpha_dual);
    double cv_curr = fatropdata_->CVL1Curr();
    double obj_curr = fatropdata_->obj_curr;
    double barrier_curr = fatropdata_->EvalBarrierCurr(mu);
    obj_curr += barrier_curr;
    double lin_decr_curr = fatropdata_->LinDecrCurr();
    double barrier_decr_curr = fatropdata_->EvalBarrierLinDecr(mu);
    lin_decr_curr += barrier_decr_curr;
    double theta_min = fatropdata_->theta_min;
    // calculation of alpha_min
    double alpha_min = gamma_alpha *
                       (lin_decr_curr > 0 ? gamma_theta
                                          : MIN(gamma_theta,
                                                cv_curr < theta_min ? MIN(-gamma_phi * cv_curr / lin_decr_curr, delta * pow(cv_curr, s_theta) / (pow(-lin_decr_curr, s_phi)))
                                                                    : -gamma_phi * cv_curr / lin_decr_curr));

    // cout << "alpha_min " << alpha_min << endl;
    // cout << "lindecr " << lin_decr_curr << endl;
    for (int ll = 1; ll < 50; ll++)
    {
        TryStep(alpha_primal, alpha_dual);
        if (alpha_primal < alpha_min)
        {
            return 0;
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
            return 1;
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
                    return ll;
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
                        (journaller_->it_curr).type = 'h';
                    }
                    fatropdata_->TakeStep();
                    journaller_->it_curr.alpha_pr = alpha_primal;
                    journaller_->it_curr.alpha_du = alpha_dual;
                    return ll;
                }
            }
        }
        alpha_primal /= 2.0;
    }
    assert(false);
    return 0;
};