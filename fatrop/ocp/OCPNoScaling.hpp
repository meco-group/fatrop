#ifndef NOSCALINGMETHODINCLUDED
#define NOSCALINGMETHODINCLUDED
#include "solver/FatropData.hpp"
#include "templates/NLPAlg.hpp"
#include "aux/SmartPtr.hpp"
#include "solver/AlgStrategy.hpp"
#include "OCPScalingMethod.hpp"
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
namespace fatrop
{
    class OCPNoScaling : public OCPScalingMethod
    {
    public:
        OCPNoScaling(const RefCountPtr<FatropParams>& fatrop_params):OCPScalingMethod(fatrop_params){};
        virtual int ComputeScalings(
            OCPKKTMemory *OCP,
            double &obj_scale,
            FatropVecBF &x_scales,
            FatropVecBF &lam_scales,
            const FatropVecBF &grad_curr) override{
                obj_scale = 1.0;
                VECSE(x_scales.nels(), 1.0, (VEC*) x_scales, 0);
                VECSE(lam_scales.nels(), 1.0, (VEC*) lam_scales, 0);
                return 0;
        };
    };

} // namespace fatrop
#endif // !SCALINGMETHODINCLUDEDLUDED