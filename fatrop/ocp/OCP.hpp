#ifndef OCPINCLUDED
#define OCPINCLUDED
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "OCPDims.hpp"
#include "aux/SmartPtr.hpp"
#include "ocp/OCPKKT.hpp"
namespace fatrop
{ 
    /** \brief interface class for OCP operations*/
    class OCP : public RefCountedObj
    {
    public:
        virtual int evalHess(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &lam) = 0;
        virtual int evalJac(
            OCPKKTMemory *OCP,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars) = 0;
        virtual int EvalConstraintViolation(
            OCPKKTMemory *OCP,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            FatropVecBF &constraint_violation) = 0;
        virtual int EvalGrad(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            FatropVecBF &gradient) = 0;
        virtual int EvalObj(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            double &res) = 0;
        virtual OCPDims GetOCPDims() const = 0;
    };
} // namespace fatrop
#endif // OCPINCLUDED