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
                                             lam_curr(nlpdims.neqs, 1),
                                             lam_next(nlpdims.neqs, 1),
                                             g_curr(nlpdims.neqs, 1),
                                             g_next(nlpdims.neqs, 1),
                                             grad_curr(nlpdims.nvars, 1),
                                             grad_next(nlpdims.nvars, 1){};
        const NLPDims nlpdims;
        FatropMemoryVecBF x_curr;
        FatropMemoryVecBF x_next;
        FatropMemoryVecBF lam_curr;
        FatropMemoryVecBF lam_next;
        FatropMemoryVecBF g_curr;
        FatropMemoryVecBF g_next;
        FatropMemoryVecBF grad_curr;
        FatropMemoryVecBF grad_next;
    };
}
#endif // FATROPDATAINCLUDED