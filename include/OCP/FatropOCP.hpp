#ifndef OCPALGINCLUDED
#define OCPALGINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "TEMPLATES/NLPAlg.hpp"
#include "AUX/SmartPtr.hpp"
#include "OCPKKT.hpp"
#include "OCPLinearSolver.hpp"
#include "AUX/FatropMemory.hpp"
#include "OCPScalingMethod.hpp"
#include "DuInfEvaluator.hpp"
#include "OCPInitializer.hpp"
#include "SPARSE/SparseOCP.hpp"
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
            const FatropVecBF &lam) override
        {
            return ocp_->evalHess(
                &ocpkktmemory_,
                obj_scale,
                primal_vars,
                lam);
        };
        int EvalJac(
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars) override
        {
            return ocp_->evalJac(
                &ocpkktmemory_,
                primal_vars,
                slack_vars);
        };
        int ComputeSD(
            const double inertia_correction_w,
            const double inertia_correction_c,
            const double mu,
            const double kappa_d,
            const FatropVecBF &dprimal_vars,
            const FatropVecBF &dlam,
            const FatropVecBF &s,
            const FatropVecBF &zL_curr,
            const FatropVecBF &zU_curr,
            const FatropVecBF &delta_zL,
            const FatropVecBF &delta_zU,
            const FatropVecBF &lower_bound,
            const FatropVecBF &upper_bound,
            const FatropVecBF &delta_s) override
        {
            // ls_ = RefCountPtr<OCPLinearSolver>(new Sparse_OCP(ocp_->GetOCPDims(), ocpkktmemory_));
            return ls_->computeSD(
                &ocpkktmemory_,
                inertia_correction_w,
                inertia_correction_c,
                mu,
                kappa_d,
                dprimal_vars,
                dlam,
                s,
                zL_curr,
                zU_curr,
                delta_zL,
                delta_zU,
                lower_bound,
                upper_bound,
                delta_s);
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
            const FatropVecBF &slack_vars,
            FatropVecBF &constraint_violation) override
        {
            return ocp_->EvalConstraintViolation(
                &ocpkktmemory_,
                primal_vars,
                slack_vars,
                constraint_violation);
        };
        int EvalGrad(
            double obj_scale,
            const FatropVecBF &primal_vars,
            FatropVecBF &gradient) override
        {

            return ocp_->EvalGrad(
                &ocpkktmemory_,
                obj_scale,
                primal_vars,
                gradient);
        };
        int EvalObj(
            double obj_scale,
            const FatropVecBF &primal_vars,
            double &res) override
        {

            return ocp_->EvalObj(
                &ocpkktmemory_,
                obj_scale,
                primal_vars,
                res);
        };
        int EvalDuInf(
            double obj_scale,
            const FatropVecBF &lam,
            const FatropVecBF &grad,
            FatropVecBF &du_inf) override
        {
            return duinfevaluator_.DuInfEval(
                &ocpkktmemory_,
                obj_scale,
                lam,
                grad,
                du_inf);
        }
        int Initialization(
            const FatropVecBF &grad,
            FatropVecBF &dlam,
            const FatropVecBF &ux_dummy,
            const FatropVecBF &s_dummy,
            FatropVecBF &s_curr,
            const FatropVecBF &zL,
            const FatropVecBF &zU,
            const FatropVecBF &lower,
            const FatropVecBF &upper) override
        {
            // assume constraint jacobian evaluated
            OCPInitializer_.AdaptKKTInitial(&ocpkktmemory_, grad, s_curr);
            return ls_->SolveInitialization(&ocpkktmemory_, dlam, ux_dummy, s_dummy, zL, zU, lower, upper);
        }

        NLPDims GetNLPDims() const override
        {
            NLPDims res;
            // states + inputs
            res.nvars = sum(ocpkktmemory_.nu) + sum(ocpkktmemory_.nx);
            // stagewise equality contraints + dynamics constraints + slack constraints
            res.neqs = sum(ocpkktmemory_.ng) + sum(ocpkktmemory_.nx) - ocpkktmemory_.nx.at(0) + sum(ocpkktmemory_.ng_ineq);
            res.nineqs = sum(ocpkktmemory_.ng_ineq);
            return res;
        };

    public:
        RefCountPtr<OCP> ocp_;
        RefCountPtr<OCPLinearSolver> ls_;
        RefCountPtr<OCPScalingMethod> scaler_;
        DuInfEvaluator duinfevaluator_;
        OCPKKTMemory ocpkktmemory_;
        OCPInitializer OCPInitializer_;
    };
} // namespace fatrop
#endif //  OCPALGINCLUDED