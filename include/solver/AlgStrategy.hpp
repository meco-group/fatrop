#ifndef ALGSTRATEGYINCLUDED
#define ALGSTRATEGYINCLUDED
#include "AUX/SmartPtr.hpp"
#include "FatropParams.hpp"
namespace fatrop
{
    class AlgStrategy: public RefCountedObj
    {
        public:
        AlgStrategy(const RefCountPtr<FatropParams>& fatrop_params):fatrop_params_(fatrop_params){};
        RefCountPtr<FatropParams> fatrop_params_;
    };
};
#endif // !ALGSTRATEGYINCLUDED