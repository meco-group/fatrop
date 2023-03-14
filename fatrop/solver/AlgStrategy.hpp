#ifndef ALGSTRATEGYINCLUDED
#define ALGSTRATEGYINCLUDED
#include "aux/SmartPtr.hpp"
#include "FatropParams.hpp"
#include <memory>
using namespace std;
namespace fatrop
{
    class AlgStrategy
    {
    public:
        AlgStrategy(const shared_ptr<FatropOptions> &fatrop_params) : fatrop_params_(fatrop_params){};
        shared_ptr<FatropOptions> fatrop_params_;
        void Initialize() {};
    };
};
#endif // !ALGSTRATEGYINCLUDED