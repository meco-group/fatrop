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
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_zL,
            const FatropVecBF &delta_zU,
            const FatropVecBF &delta_s,
            const FatropVecBF &sigma_L,
            const FatropVecBF &sigma_U,
            const FatropVecBF &gradb_L,
            const FatropVecBF &gradb_U,
            const FatropVecBF &gradb_plus,
            const FatropVecBF &zL_curr,
            const FatropVecBF &zU_curr) override;
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
            FatropVecBF &s_curr,
            const FatropVecBF &zL,
            const FatropVecBF &zU) override;

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
        FatropMemoryVecBF s_memvec;
        FatropMemoryVecBF ux_memvec;
        FatropVecBF sigma;
        FatropVecBF gradb;
        FatropVecBF s_dummy;
        FatropVecBF s_zero;
        FatropVecBF ux_dummy;
    };
} // namespace fatrop
#endif //  OCPALGINCLUDED