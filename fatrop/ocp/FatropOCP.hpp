#ifndef OCPALGINCLUDED
#define OCPALGINCLUDED
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "templates/NLPAlg.hpp"
#include "aux/SmartPtr.hpp"
#include "OCPKKT.hpp"
#include "OCPLinearSolver.hpp"
#include "OCPScalingMethod.hpp"
#include "DuInfEvaluator.hpp"
#include "OCPInitializer.hpp"
#include "ocp/LineSearchDDP.hpp"
// #include "sparse/SparseOCP.hpp"
#include "OCP.hpp"
#include <memory>
using namespace std;
// #include <unistd.h>
namespace fatrop
{
    class FatropOCP : public FatropNLP
    {
    public:
        FatropOCP(
            const shared_ptr<OCP> &ocp,
            const shared_ptr<OCPLinearSolver> &ls,
            const shared_ptr<OCPScalingMethod> &scaler);
        int EvalHess(
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &lam) override;
        int EvalJac(
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars) override;
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
            const FatropVecBF &delta_s) override;
        int ComputeScalings(
            double &obj_scale,
            FatropVecBF &x_scales,
            FatropVecBF &lam_scales,
            const FatropVecBF &grad_curr);
        int EvalConstraintViolation(
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            FatropVecBF &constraint_violation) override;
        int EvalGrad(
            double obj_scale,
            const FatropVecBF &primal_vars,
            FatropVecBF &gradient) override;
        int EvalObj(
            double obj_scale,
            const FatropVecBF &primal_vars,
            double &res) override;
        int EvalDuInf(
            double obj_scale,
            const FatropVecBF &lam,
            const FatropVecBF &grad,
            FatropVecBF &du_inf) override;
        int Initialization(
            const FatropVecBF &grad,
            FatropVecBF &dlam,
            const FatropVecBF &ux_dummy,
            const FatropVecBF &s_dummy,
            FatropVecBF &s_curr,
            const FatropVecBF &zL,
            const FatropVecBF &zU,
            const FatropVecBF &lower,
            const FatropVecBF &upper) override;

        NLPDims GetNLPDims() const override;
        void Finalize() override;
        void Reset() override;

    public:
        shared_ptr<OCP> ocp_;
        shared_ptr<OCPLinearSolver> ls_;
        shared_ptr<OCPScalingMethod> scaler_;
        DuInfEvaluator duinfevaluator_;
        OCPKKTMemory ocpkktmemory_;
        OCPInitializer OCPInitializer_;
        double FE_time = 0.0;
        double hess_time = 0.0;
        double jac_time = 0.0;
        double cv_time = 0.0;
        double grad_time = 0.0;
        double obj_time = 0.0;
        int hess_count = 0;
        int jac_count = 0;
        int cv_count = 0;
        int grad_count = 0;
        int obj_count = 0;
    };
} // namespace fatrop
#endif //  OCPALGINCLUDED