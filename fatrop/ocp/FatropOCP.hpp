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
#include "solver/FatropPrinter.hpp"
// #include "ocp/LineSearchDDP.hpp"
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
        // FatropOCP();
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
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_s,
            const FatropVecBF &sigma_total,
            const FatropVecBF &gradb_total) override;
        int SolveSOC(
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_s,
            const FatropVecBF &constraint_violation) override;
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
        int Initialization_s(
            FatropVecBF &s_curr) override;
        int Initialization_dual(
            const FatropVecBF &grad,
            FatropVecBF &dlam,
            const FatropVecBF &zL,
            const FatropVecBF &zU) override;

        int GetBounds(
            FatropVecBF &lower,
            FatropVecBF &upper) const override
        {
            return ocp_->GetBounds(lower, upper);
        };
        int GetInitialGuess(
            FatropVecBF &initial) const override
        {
            return ocp_->GetInitialGuess(initial);
        };
        // int GetDefaultParams(
        //     FatropParams &params) const override
        //     {
        //        return ocp_->GetDefaultParams(params);
        //     };
        NLPDims GetNLPDims() const override;
        void Finalize() override;
        void Reset() override;

    public:
        shared_ptr<OCP> ocp_;
        OCPDims dims_;
        NLPDims nlpdims_;
        shared_ptr<OCPLinearSolver> ls_;
        shared_ptr<OCPScalingMethod> scaler_;
        DuInfEvaluator duinfevaluator_;
        OCPKKTMemory ocpkktmemory_;
        OCPInitializer OCPInitializer_;
        FatropMemoryVecBF s_memvec;
        FatropMemoryVecBF ux_memvec;
        FatropVecBF sigma;
        FatropVecBF gradb;
        FatropVecBF s_dummy;
        FatropVecBF s_zero;
        FatropVecBF ux_dummy;
        FatropMemoryVecBF rhs_rq;
        FatropMemoryVecBF rhs_b;
        FatropMemoryVecBF rhs_g;
        FatropMemoryVecBF rhs_g_ineq;
        FatropMemoryVecBF rhs_gradb;
        FatropMemoryVecBF rhs_rq2;
        FatropMemoryVecBF rhs_b2;
        FatropMemoryVecBF rhs_g2;
        FatropMemoryVecBF rhs_g_ineq2;
        FatropMemoryVecBF rhs_gradb2;
        FatropMemoryVecBF gradb_total_cache;
        FatropMemoryVecBF sigma_total_cache;
        FatropMemoryVecBF ux_test;
        FatropMemoryVecBF lam_test;
        FatropMemoryVecBF delta_s_test;
        double inertia_correction_w_cache;
        double inertia_correction_c_cache;
    };
} // namespace fatrop
#endif //  OCPALGINCLUDED