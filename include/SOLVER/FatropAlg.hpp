#ifndef FATROPALGINCLUDED
#define FATROPALGINCLUDED
#include "AUX/SmartPtr.hpp"
#include "TEMPLATES/NLPAlg.hpp"
#include "FatropData.hpp"
#include "Filter.hpp"
#include "LineSearch.hpp"
namespace fatrop
{
    struct IterationData
    {
        
    };
    class FatropAlg: public RefCountedObj
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
            return 0.0;
        }
        double EvalObjNext()
        {
            return 0.0;
        }
        int ComputeSD()
        {
            return 0;
        }
        RefCountPtr<FatropNLP> fatropnlp_;
        RefCountPtr<FatropData> fatropdata_;
    };
} // namespace fatrop
#endif // FATROPALGINCLUDED