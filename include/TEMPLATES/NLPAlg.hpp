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
    class FatropNLP : public RefCountedObj
    {
    public:
        virtual int EvalHess(
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &lam) = 0;
        virtual int EvalJac(
            const FatropVecBF &primal_vars) = 0;
        virtual int EvalConstraintViolation(
            const FatropVecBF &primal_vars,
            FatropVecBF &constraint_violation) = 0;
        virtual int EvalGrad(
            double obj_scale,
            const FatropVecBF &primal_vars,
            FatropVecBF &gradient) = 0;
        virtual int EvalObj(
            double obj_scale,
            const FatropVecBF &primal_vars,
            double &res) = 0;
        virtual int EvalDuInf(
            double obj_scale,
            const FatropVecBF &lam,
            const FatropVecBF &grad,
            FatropVecBF &du_inf) = 0;
        virtual int ComputeSD(
            const double intertia_correction,
            const FatropVecBF &dprimal_vars,
            const FatropVecBF &dlam) = 0;
        virtual NLPDims GetNLPDims() const = 0;
        virtual int ComputeScalings(
            double &obj_scale,
            FatropVecBF &x_scales,
            FatropVecBF &lam_scales,
            const FatropVecBF &grad_curr) = 0;
        virtual int Initializiaton(
            const FatropVecBF& grad,
            FatropVecBF& lam,
            FatropVecBF& optimvarsdummy) = 0 ;
    };
} // namespace fatrop
#endif // NLPINCLUDED