// solver data
#ifndef FATROPDATAINCLUDED
#define FATROPDATAINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "TEMPLATES/NLPAlg.hpp"
using namespace std;
namespace fatrop{
struct FatropData: public RefCountedObj
{
    const NLPDims nlpdims;
    FatropVecBF x_curr;
    FatropVecBF x_next;
    FatropVecBF lam_curr;
    FatropVecBF lam_next;
    FatropVecBF g_curr;
    FatropVecBF g_next;
    FatropVecBF grad_curr;
    FatropVecBF grad_next;
};
}
#endif // FATROPDATAINCLUDED