#ifndef SCALINGMETHODINCLUDED
#define SCALINGMETHODINCLUDED
#include "SOLVER/FatropData.hpp"
#include "TEMPLATES/NLPAlg.hpp"
#include "AUX/SmartPtr.hpp"
#include "SOLVER/AlgStrategy.hpp"
namespace fatrop
{
    class OCPScalingMethod  //this is an OCP strategy
    {
    public:
        virtual int ComputeScalings(
            OCPKKTMemory *OCP,
            double& obj_scale,
            FatropVecBF &x_scales,
            FatropVecBF &lam_scales,
            const FatropVecBF &grad_curr) = 0;
    };

} // namespace fatrop
#endif // !SCALINGMETHODINCLUDED