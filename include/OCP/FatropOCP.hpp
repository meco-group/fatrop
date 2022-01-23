#ifndef OCPALGINCLUDED
#define OCPALGINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "TEMPLATES/NLPAlg.hpp"
#include "AUX/SmartPtr.hpp"
#include "OCPKKT.hpp"
#include "OCPLinearSolver.hpp"
#include "AUX/FatropMemory.hpp"
#include "OCPScalingMethod.hpp"
namespace fatrop
{
    class FatropOCP : public FatropNLP
    {
    public:
        FatropOCP(
            const RefCountPtr<OCP> &ocp,
            const RefCountPtr<OCPLinearSolver> &ls,
            const RefCountPtr<OCPScalingMethod> &scaler) : ocp_(ocp), ls_(ls), scaler_(scaler), ocpkktmemory_(ocp_->GetOCPDims()){};
        int EvalHess(
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &scales_primal_vars,
            const FatropVecBF &lam,
            const FatropVecBF &scales_lam) override
        {
            return ocp_->evalHess(
                &ocpkktmemory_,
                obj_scale,
                primal_vars,
                scales_primal_vars,
                lam,
                scales_lam);
        };
        int EvalJac(
            const FatropVecBF &primal_vars,
            const FatropVecBF &scales_primal_vars,
            const FatropVecBF &scales_lam) override
        {
            return ocp_->evalJac(
                &ocpkktmemory_,
                primal_vars,
                scales_primal_vars,
                scales_lam);
        };
        int ComputeSD(
            const double inertia_correction,
            const FatropVecBF &dprimal_vars,
            const FatropVecBF &dlam) override
        {
            return ls_->computeSD(
                &ocpkktmemory_,
                inertia_correction,
                dprimal_vars,
                dlam);
        };
        int ComputeScalings(
            double &obj_scale,
            FatropVecBF &x_scales,
            FatropVecBF &lam_scales,
            const FatropVecBF &grad_curr)
        {
            return scaler_->ComputeScalings(
                &ocpkktmemory_,
                obj_scale,
                x_scales,
                lam_scales,
                grad_curr);
        };
        int EvalConstraintViolation(
            const FatropVecBF &primal_vars,
            const FatropVecBF &scales_primal_vars,
            const FatropVecBF &scales_lam,
            const FatropVecBF &constraint_violation) override
        {
            return ocp_->EvalConstraintViolation(
                primal_vars,
                scales_primal_vars,
                scales_lam,
                constraint_violation);
        };
        NLPDims GetNLPDims() const override
        {
            NLPDims res;
            // states + inputs
            res.nvars = sum(ocpkktmemory_.nu) + sum(ocpkktmemory_.nx);
            // stagewise equality contraints + dynamics constraints
            res.neqs = sum(ocpkktmemory_.ng) + sum(ocpkktmemory_.nx) - ocpkktmemory_.nx.at(0);
            return res;
        };

    public:
        RefCountPtr<OCP> ocp_;
        RefCountPtr<OCPLinearSolver> ls_;
        RefCountPtr<OCPScalingMethod> scaler_;
        OCPKKTMemory ocpkktmemory_;
    };
} // namespace fatrop
#endif //  OCPALGINCLUDED