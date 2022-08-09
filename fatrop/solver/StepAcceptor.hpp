#ifndef STEPACCEPTORINCLUDED
#define STEPACCEPTORINCLUDED
#include "aux/SmartPtr.hpp"
#include "FatropData.hpp"
#include "FatropParams.hpp"
#include "Filter.hpp"
#include <memory>
using namespace std;
namespace fatrop
{
    class StepAcceptor
    {
        StepAcceptor(
        const shared_ptr<FatropData>& fatropdata,
        const shared_ptr<Filter>& filter,
        const shared_ptr<FatropParams>& fatropparams
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
        shared_ptr<FatropData> fatropdata_;
        shared_ptr<Filter> filter_;
        shared_ptr<FatropParams> fatropparams_;
    };
} // namespace fatrop

#endif // STEPACCEPTORINCLUDED