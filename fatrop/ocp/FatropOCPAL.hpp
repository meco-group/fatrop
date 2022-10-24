#ifndef FATROPOCPALINCLUDED
#define FATROPOCPALINCLUDED
#include "FatropOCP.hpp"
#include "templates/FatropNLPAL.hpp"
#include "ocp/OCPAL.hpp"
using namespace std;
namespace fatrop
{
    class FatropOCPAL : public FatropOCP, public FatropNLPAL
    {
        public:
        FatropOCPAL(
            const shared_ptr<OCPAL> &ocp,
            const shared_ptr<OCPLinearSolver> &ls,
            const shared_ptr<OCPScalingMethod> &scaler) : FatropOCP(ocp, ls, scaler), ocp_(ocp)
        {
        }
        int EvalHess(
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &lam) override
        {
            return FatropOCP::EvalHess(
                obj_scale,
                primal_vars,
                lam);
        }
        int EvalJac(
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars) override
        {
            return FatropOCP::EvalJac(
                primal_vars,
                slack_vars);
        }
        int ComputeSD(
            const double inertia_correction_w,
            const double inertia_correction_c,
            const double mu,
            const double kappa_d,
            const FatropVecBF &dprimal_vars,
            const FatropVecBF &dlam,
            const FatropVecBF &lam_curr,
            const FatropVecBF &s,
            const FatropVecBF &zL_curr,
            const FatropVecBF &zU_curr,
            const FatropVecBF &delta_zL,
            const FatropVecBF &delta_zU,
            const FatropVecBF &lower_bound,
            const FatropVecBF &upper_bound,
            const FatropVecBF &delta_s) override
        {
            return FatropOCP::ComputeSD(
                inertia_correction_w,
                inertia_correction_c,
                mu,
                kappa_d,
                dprimal_vars,
                dlam,
                lam_curr,
                s,
                zL_curr,
                zU_curr,
                delta_zL,
                delta_zU,
                lower_bound,
                upper_bound,
                delta_s);
        }
        int ComputeResidual(
            const double inertia_correction_w,
            const double inertia_correction_c,
            const double mu,
            const double kappa_d,
            const FatropVecBF &dprimal_vars,
            const FatropVecBF &dlam,
            const FatropVecBF &lam_curr,
            const FatropVecBF &s,
            const FatropVecBF &zL_curr,
            const FatropVecBF &zU_curr,
            const FatropVecBF &delta_zL,
            const FatropVecBF &delta_zU,
            const FatropVecBF &lower_bound,
            const FatropVecBF &upper_bound,
            const FatropVecBF &delta_s,
            FatropVecBF &residual) override
        {
            return FatropOCP::ComputeResidual(
                inertia_correction_w,
                inertia_correction_c,
                mu,
                kappa_d,
                dprimal_vars,
                dlam,
                lam_curr,
                s,
                zL_curr,
                zU_curr,
                delta_zL,
                delta_zU,
                lower_bound,
                upper_bound,
                delta_s, 
                residual);
        }
        int ComputeScalings(
            double &obj_scale,
            FatropVecBF &x_scales,
            FatropVecBF &lam_scales,
            const FatropVecBF &grad_curr) override
        {
            return FatropOCP::ComputeScalings(
                obj_scale,
                x_scales,
                lam_scales,
                grad_curr);
        }
        int EvalConstraintViolation(
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            FatropVecBF &constraint_violation) override
        {
            return FatropOCP::EvalConstraintViolation(
                primal_vars,
                slack_vars,
                constraint_violation);
        }
        int EvalGrad(
            double obj_scale,
            const FatropVecBF &primal_vars,
            FatropVecBF &gradient) override
        {
            return FatropOCP::EvalGrad(
                obj_scale,
                primal_vars,
                gradient);
        }
        int EvalObj(
            double obj_scale,
            const FatropVecBF &primal_vars,
            double &res) override
        {
            return FatropOCP::EvalObj(
                obj_scale,
                primal_vars,
                res);
        }
        int EvalDuInf(
            double obj_scale,
            const FatropVecBF &lam,
            const FatropVecBF &grad,
            FatropVecBF &du_inf) override
        {
            return FatropOCP::EvalDuInf(
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
            return FatropOCP::Initialization(
                grad,
                dlam,
                ux_dummy,
                s_dummy,
                s_curr,
                zL,
                zU,
                lower,
                upper);
        }

        NLPDims GetNLPDims() const override
        {
            return FatropOCP::GetNLPDims();
        }

        void Reset() override
        {
            return FatropOCP::Reset();
        }
        int SetIneqsBounds(const vector<double> &lower_boundsin, const vector<double> &upper_boundsin)
        {
            return ocp_->SetIneqsBounds(lower_boundsin, upper_boundsin);
        };
        int EvalInequalities(
            const FatropVecBF &primal_vars,
            FatropVecBF &inequalities)
        {
            return ocp_->EvalInequalities(
                &ocpkktmemory_,
                primal_vars,
                inequalities);
        }
        int SetIneqLagrMult(const FatropVecBF &ineqlagrmultL, const FatropVecBF &ineqlagrmultU)
        {
            return ocp_->SetIneqLagrMult(ineqlagrmultL, ineqlagrmultU);
        }
        int ResetIneqLagrMult()
        {
            return ocp_->ResetIneqLagrMult();
        }
        int SetPenalty(double penalty)
        {
            return ocp_->SetPenalty(penalty);
        };

        int GetNOIneqs()
        {
            return ocp_->GetTotalNOIneqs();
        };
        const shared_ptr<OCPAL> ocp_;
    };

} // namespace fatrop
#endif // FATROPOCPALINCLUDED