#ifndef SCALINGMETHODINCLUDED
#define SCALINGMETHODINCLUDED
#include "FatropData.hpp"
#include "TEMPLATES/NLPAlg.hpp"
#include "AUX/SmartPtr.hpp"
#include "AlgStrategy.hpp"
namespace fatrop
{
    class ScalingMethod : public AlgStrategy
    {
        ScalingMethod(const RefCountPtr<FatropData> &fatropdata) : fatropdata_(fatropdata){};
        RefCountPtr<FatropData> fatropdata_;
    };
} // namespace fatrop
#endif // !SCALINGMETHODINCLUDED