#ifndef FATROPALGINCLUDED
#define FATROPALGINCLUDED
#include "AUX/SmartPtr.hpp"
#include "TEMPLATES/NLPAlg.hpp"
#include "FatropData.hpp"
#include "Filter.hpp"
namespace fatrop
{
    class FatropAlg
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
            // fatropnlp_->EvalHess();
            // fatropnlp_->EvalJac();
            // fatropnlp_->ComputeSD();
            return 0;
        }
        int EvalHess()
        {
            // return fatropnlp_->EvalHess
            // (
            // ) 
            return 0;
        }
        int EvalJac()
        {
            return 0;
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