#ifndef NLPINCLUDED
#define NLPINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
namespace fatrop
{
    struct NLPDims
    {
        int nvars;
        int neqs;
    };
    class NLPAlg
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
    };
} // namespace fatrop
#endif // NLPINCLUDED