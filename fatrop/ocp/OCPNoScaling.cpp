#include "ocp/OCPNoScaling.hpp"
using namespace fatrop;
OCPNoScaling::OCPNoScaling(const shared_ptr<FatropParams> &fatrop_params) : OCPScalingMethod(fatrop_params){};


int OCPNoScaling::ComputeScalings(
    OCPKKTMemory *OCP,
    double &obj_scale,
    FatropVecBF &x_scales,
    FatropVecBF &lam_scales,
    const FatropVecBF &grad_curr) 
{
    obj_scale = 1.0;
    VECSE(x_scales.nels(), 1.0, (VEC *)x_scales, 0);
    VECSE(lam_scales.nels(), 1.0, (VEC *)lam_scales, 0);
    return 0;
};