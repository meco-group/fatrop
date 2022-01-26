#ifndef FATROPALGINCLUDED
#define FATROPALGINCLUDED
#include "AUX/SmartPtr.hpp"
#include "TEMPLATES/NLPAlg.hpp"
#include "FatropData.hpp"
#include "Filter.hpp"
#include "LineSearch.hpp"
#include "StepAcceptor.hpp"
namespace fatrop
{
    struct IterationData
    {
    };
    class FatropAlg : public RefCountedObj
    {
    public:
        FatropAlg(
            const RefCountPtr<FatropNLP> &fatropnlp,
            const RefCountPtr<FatropData> &fatropdata) : fatropnlp_(fatropnlp),
                                                         fatropdata_(fatropdata)
        {
        }
        int Optimize()
        {
            return 0;
        }
        inline int EvalHess()
        {
            return fatropnlp_->EvalHess(
                fatropdata_->obj_scale,
                fatropdata_->x_curr,
                fatropdata_->lam_curr);
        }
        inline int EvalJac()
        {
            return fatropnlp_->EvalJac(
                fatropdata_->x_curr);
        }
        inline int EvalCVCurr()
        {
            return fatropnlp_->EvalConstraintViolation(
                fatropdata_->x_curr,
                fatropdata_->g_curr);
        }
        inline int EvalCVNext()
        {
            return fatropnlp_->EvalConstraintViolation(
                fatropdata_->x_next,
                fatropdata_->g_next);
        }
        inline int EvalGradCurr()
        {
            return fatropnlp_->EvalGrad(
                fatropdata_->obj_scale,
                fatropdata_->x_curr,
                fatropdata_->grad_curr);
        }
        int EvalGradNext()
        {
            return fatropnlp_->EvalGrad(
                fatropdata_->obj_scale,
                fatropdata_->x_next,
                fatropdata_->grad_next);
        }
        double EvalObjCurr()
        {
            double res = 0.0;
            fatropnlp_->EvalObj(
                fatropdata_->obj_scale,
                fatropdata_ -> x_curr,
                res);
            return res;
        }
        double EvalObjNext()
        {
            double res = 0.0;
            fatropnlp_->EvalObj(
                fatropdata_->obj_scale,
                fatropdata_ -> x_next,
                res);
            return res;
        }
        int ComputeSD(double inertia_correction)
        {
            return fatropnlp_->ComputeSD(
                inertia_correction,
                fatropdata_->delta_x,
                fatropdata_->lam_calc);
        }
        RefCountPtr<FatropNLP> fatropnlp_;
        RefCountPtr<FatropData> fatropdata_;
    };
} // namespace fatrop
#endif // FATROPALGINCLUDED