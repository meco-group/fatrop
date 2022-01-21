// solver data
#ifndef FATROPDATAINCLUDED
#define FATROPDATAINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "TEMPLATES/NLPAlg.hpp"
using namespace std;
namespace fatrop
{
    struct FatropData : public RefCountedObj
    {
        FatropData(const NLPDims &nlpdims) : nlpdims(nlpdims),
                                             x_curr(nlpdims.nvars, 1),
                                             x_next(nlpdims.nvars, 1),
                                             delta_x(nlpdims.nvars, 1),
                                             x_scales(nlpdims.nvars, 1),
                                             lam_curr(nlpdims.neqs, 1),
                                             lam_next(nlpdims.neqs, 1),
                                             lam_scales(nlpdims.neqs, 1),
                                             g_curr(nlpdims.neqs, 1),
                                             g_next(nlpdims.neqs, 1),
                                             grad_curr(nlpdims.nvars, 1),
                                             grad_next(nlpdims.nvars, 1){};
        const NLPDims nlpdims;
        double obj_scaling_factor = 1.0;
        FatropMemoryVecBF x_curr;
        FatropMemoryVecBF x_next;
        FatropMemoryVecBF delta_x;
        FatropMemoryVecBF x_scales;
        FatropMemoryVecBF lam_curr;
        FatropMemoryVecBF lam_next;
        FatropMemoryVecBF lam_scales;
        FatropMemoryVecBF g_curr;
        FatropMemoryVecBF g_next;
        FatropMemoryVecBF grad_curr;
        FatropMemoryVecBF grad_next;
    };
}
#endif // FATROPDATAINCLUDED