#ifndef FATROPNLPALINCLUDED
#define FATROPNLPALINCLUDED
#include "NLPAlg.hpp"
namespace fatrop
{
    class FatropNLPAL : public FatropNLP
    {
    public:
        virtual int EvalInequalities(
            const FatropVecBF &primal_vars,
            FatropVecBF &inequalities) = 0;
        virtual int SetIneqsBounds(const vector<double> &lower_boundsin, const vector<double> &upper_boundsin) = 0;
        virtual int SetIneqLagrMult(const FatropVecBF &ineqlagrmult) = 0;
        virtual int SetPenalty(double penalty) = 0;
    };
} // namespace fatrop
#endif // FATROPNLPALINCLUDED