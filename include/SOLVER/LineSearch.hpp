#ifndef LINESEARCHINCLUDED
#define LINESEARCHINCLUDED
#include "AlgStrategy.hpp"
#include "FatropAlg.hpp"
#include <cmath>
namespace fatrop
{
    class LineSearch : public AlgStrategy
    {
    public:
        LineSearch(
            const RefCountPtr<FatropParams> &fatropparams,
            const RefCountPtr<FatropNLP> &nlp,
            const RefCountPtr<FatropData> &fatropdata) : AlgStrategy(fatropparams),
                                                         fatropnlp_(nlp),
                                                         fatropdata_(fatropdata){};
        virtual int FindAcceptableTrialPoint()=0;
        inline int EvalCVNext()
        {
            return fatropnlp_->EvalConstraintViolation(
                fatropdata_->x_next,
                fatropdata_->g_next);
        }
        double EvalObjNext()
        {
            double res = 0.0;
            fatropnlp_->EvalObj(
                fatropdata_->obj_scale,
                fatropdata_->x_next,
                res);
            return res;
        }
        RefCountPtr<FatropNLP> fatropnlp_;
        RefCountPtr<FatropData> fatropdata_;
    };

    class BackTrackingLineSearch : public LineSearch
    {
    public:
        BackTrackingLineSearch(
            const RefCountPtr<FatropParams> &fatropparams,
            const RefCountPtr<FatropNLP> &nlp,
            const RefCountPtr<FatropData> &fatropdata,
            const RefCountPtr<Filter> &filter)
            : LineSearch(fatropparams, nlp, fatropdata), filter_(filter){
                Initialize();
            };
        void Initialize()
        {
            // todo avoid reallocation when maxiter doesn't change
            s_phi = fatrop_params_->s_phi;
            delta = fatrop_params_->delta;
            s_theta = fatrop_params_->s_theta;
            gamma_theta = fatrop_params_->gamma_theta;
            gamma_phi = fatrop_params_->gamma_phi;
            eta_phi = fatrop_params_->eta_phi;
        }
        int FindAcceptableTrialPoint()
        {
            double alpha_primal = 1.0;
            double alpha_dual = 1.0;
            double cv_curr = fatropdata_->CVL1Curr();
            double obj_curr = fatropdata_->obj_curr;
            double lin_decr_curr = fatropdata_->LinDecrCurr();
            for (int ll = 0; ll < 50; ll++)
            {
                fatropdata_->TryStep(alpha_primal, alpha_dual);
                double obj_next = EvalObjNext();
                EvalCVNext();
                double cv_next = fatropdata_->CVL1Next();
                // todo change iteration number from zero to real iteration number
                if (filter_->IsAcceptable(FilterData(0, obj_next, cv_next)))
                {
                    bool switch_cond = (lin_decr_curr < 0) && (alpha_primal * pow(-lin_decr_curr, s_phi) > delta * pow(cv_curr, s_theta));
                    if (switch_cond && (cv_curr <= fatropdata_->theta_min))
                    {
                        // f-step
                        if (obj_next - obj_curr < eta_phi * alpha_primal * lin_decr_curr)
                        {
                            fatropdata_->TakeStep();
                            return ll;
                        }
                    }
                    else
                    {
                        // h-step
                        // check sufficient decrease wrt current iterate
                        if ((cv_next < (1.0 - gamma_theta) * cv_curr) || (obj_next < obj_curr - gamma_phi * cv_curr))
                        {
                            if (!switch_cond || !(obj_next - obj_curr < eta_phi * alpha_primal * lin_decr_curr))
                            {
                                filter_->Augment(FilterData(0, obj_curr - gamma_phi * cv_curr, cv_curr * (1 - gamma_theta)));
                            }
                            fatropdata_->TakeStep();
                            return ll;
                        }
                    }
                }
                alpha_primal /= 2.0;
            }
            assert(false);
            return 0;
        };
        RefCountPtr<Filter> filter_;
        double s_phi;
        double delta;
        double s_theta;
        double gamma_theta;
        double gamma_phi;
        double eta_phi;
    };
} // namespace fatrop
#endif // LINESEARCHINCLUDED