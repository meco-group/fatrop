#ifndef FATROPNLPALINCLUDED
#define FATROPNLPALINCLUDED
#include "NLPAlg.hpp"
namespace fatrop
{
    class FatropNLPAL : public FatropNLP
    {
        virtual int EvalInequalities(
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            FatropVecBF &inequalities) = 0;
    };
} // namespace fatrop
#endif // FATROPNLPALINCLUDED