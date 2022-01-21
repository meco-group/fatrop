#ifndef NOSCALINGMETHODINCLUDED
#define NOSCALINGMETHODINCLUDED
#include "SOLVER/FatropData.hpp"
#include "TEMPLATES/NLPAlg.hpp"
#include "AUX/SmartPtr.hpp"
#include "SOLVER/AlgStrategy.hpp"
#include "OCPScalingMethod.hpp"
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
namespace fatrop
{
    class OCPNoScaling : public OCPScalingMethod
    {
    public:
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