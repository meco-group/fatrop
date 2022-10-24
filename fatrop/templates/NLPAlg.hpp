#ifndef NLPINCLUDED
#define NLPINCLUDED
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "aux/SmartPtr.hpp"
namespace fatrop
{
    struct NLPDims
    {
        int nvars;
        int neqs;
        int nineqs;
    };
    class FatropNLP 
    {
    public:
        virtual int EvalHess(
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &lam) = 0;
        virtual int EvalJac(
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars) = 0;
        virtual int EvalConstraintViolation(
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
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
            const double intertia_correction_w,
            const double intertia_correction_c,
            const double mu,
            const double kappa_d,
            const FatropVecBF &dprimal_vars,
            const FatropVecBF &dlam,
            const FatropVecBF &lam_curr,
            const FatropVecBF &s,
            const FatropVecBF &zL_curr,
            const FatropVecBF &zU_curr,
            const FatropVecBF &delta_zL,
            const FatropVecBF &delta_ZU,
            const FatropVecBF &lower_bound,
            const FatropVecBF &upper_bound,
            const FatropVecBF &delta_s) = 0;
        virtual NLPDims GetNLPDims() const = 0;
        virtual int ComputeScalings(
            double &obj_scale,
            FatropVecBF &x_scales,
            FatropVecBF &lam_scales,
            const FatropVecBF &grad_curr) = 0;
        virtual int Initialization(
            const FatropVecBF &grad,
            FatropVecBF &dlam,
            const FatropVecBF &ux_dummy,
            const FatropVecBF &s_dummy,
            FatropVecBF &s_curr,
            const FatropVecBF &zL,
            const FatropVecBF &zU,
            const FatropVecBF &lower,
            const FatropVecBF &upper) = 0;
        virtual void Finalize() {};
        virtual void Reset() {};
    };
} // namespace fatrop
#endif // NLPINCLUDED