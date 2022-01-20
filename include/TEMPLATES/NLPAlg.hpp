#ifndef NLPINCLUDED
#define NLPINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "AUX/SmartPtr.hpp"
namespace fatrop
{
    struct NLPDims
    {
        int nvars;
        int neqs;
    };
    class FatropNLP: public RefCountedObj
    {
    public:
        virtual int EvalHess(
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &scales_primal_vars,
            const FatropVecBF &lam,
            const FatropVecBF &scales_lam) = 0;
        virtual int EvalJac(
            const FatropVecBF &primal_vars,
            const FatropVecBF &scales_primal_vars,
            const FatropVecBF &scales_lam) = 0;
        virtual int ComputeSD(
            const double intertia_correction,
            const FatropVecBF &dprimal_vars,
            const FatropVecBF &dlam) = 0;
            virtual NLPDims GetNLPDims() const =0;
    };
} // namespace fatrop
#endif // NLPINCLUDED