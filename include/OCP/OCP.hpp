#ifndef OCPINCLUDED
#define OCPINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "OCPDims.hpp"
#include "AUX/SmartPtr.hpp"
namespace fatrop
{
    class OCP: public RefCountedObj
    {
        public:
        virtual int evalHess(
            OCPKKTMemory* OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &scales_primal_vars,
            const FatropVecBF &lam,
            const FatropVecBF &scales_lam) = 0;
        virtual int evalJac(
            OCPKKTMemory* OCP,
            const FatropVecBF &primal_vars,
            const FatropVecBF &scales_primal_vars,
            const FatropVecBF &scales_lam) = 0;
        virtual int EvalConstraintViolation(
            OCPKKTMemory *OCP,
            const FatropVecBF &primal_vars,
            const FatropVecBF &scales_primal_vars,
            const FatropVecBF &scales_lam,
            const FatropVecBF &constraint_violation) = 0;
        virtual int EvalGrad(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &scales_primal_vars,
            const FatropVecBF &gradient) = 0;
        virtual OCPDims GetOCPDims() const  = 0;
    };
} // namespace fatrop
#endif // OCPINCLUDED