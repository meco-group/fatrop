#ifndef NLPINCLUDED
#define NLPINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
namespace fatrop
{
    class NLP
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
    };
} // namespace fatrop
#endif // NLPINCLUDED