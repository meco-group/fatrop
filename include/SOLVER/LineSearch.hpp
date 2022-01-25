#ifndef LINESEARCHINCLUDED
#define LINESEARCHINCLUDED
#include "AlgStrategy.hpp"
#include "FatropAlg.hpp"
namespace fatrop
{
    class LineSearch : public AlgStrategy
    {
    public:
        LineSearch(
            const RefCountPtr<FatropParams> &fatropparams,
            const RefCountPtr<FatropParams> &nlp,
            const RefCountPtr<FatropData> &fatropdata) : AlgStrategy(fatropparams),
                                                         nlp(nlp),
                                                         fatropdata(fatropdata){};
        RefCountPtr<FatropParams> nlp;
        RefCountPtr<FatropData> fatropdata;
    };

    class BackTrackingLineSearch : public LineSearch
    {
    public:
        BackTrackingLineSearch(
            const RefCountPtr<FatropParams> &fatropparams,
            const RefCountPtr<FatropParams> &nlp,
            const RefCountPtr<FatropData> &FatropData) : LineSearch(fatropparams, nlp, fatropdata){};
        int FindAcceptableTrialPoint(){
            return 0;
        };
    };
} // namespace fatrop
#endif // LINESEARCHINCLUDED