#ifndef DUINFEVALINCLUDED
#define DUINFEVALINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "OCPKKT.hpp"
namespace fatrop
{
    class DuInfEvaluator
    {
        public:
        int DuInfEval(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &lam)
        {
            return 0;
        }
    };
}
#endif //  DUINFEVALINCLUDED