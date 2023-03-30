#ifndef DUINFEVALINCLUDED
#define DUINFEVALINCLUDED
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "OCPKKT.hpp"
#include "aux/Common.hpp"
namespace fatrop
{
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP->name)
    class DuInfEvaluator
    {
    public:
        int evaluate(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &lam,
            const FatropVecBF &grad,
            FatropVecBF &du_inf);
    };
}
#endif //  DUINFEVALINCLUDED