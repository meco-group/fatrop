#ifndef STEPACCEPTORINCLUDED
#define STEPACCEPTORINCLUDED
#include "AUX/SmartPtr.hpp"
#include "FatropData.hpp"
#include "FatropParams.hpp"
#include "Filter.hpp"
namespace fatrop
{
    class StepAcceptor
    {
        StepAcceptor(
        const RefCountPtr<FatropData>& fatropdata,
        const RefCountPtr<Filter>& filter,
        const RefCountPtr<FatropParams>& fatropparams
    ): fatropdata_(fatropdata), filter_(filter), fatropparams_(fatropparams)
    {
        Initialize();
    };
        public:
        void Initialize()
        {

        }
        void AcceptTrialStep()
        {
        }
        private:
        RefCountPtr<FatropData> fatropdata_;
        RefCountPtr<Filter> filter_;
        RefCountPtr<FatropParams> fatropparams_;
    };
} // namespace fatrop

#endif // STEPACCEPTORINCLUDED