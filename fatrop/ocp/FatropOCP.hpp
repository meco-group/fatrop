#ifndef OCPALGINCLUDED
#define OCPALGINCLUDED
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "templates/NLPAlg.hpp"
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
// #include <unistd.h>
namespace fatrop
{
    class FatropOCP : public FatropNLP
    {
    public:
        // FatropOCP();
        FatropOCP(
            const std::shared_ptr<OCP> &ocp,
            const std::shared_ptr<OCPLinearSolver> &ls,
            const std::shared_ptr<OCPScalingMethod> &scaler, const std::shared_ptr<FatropOptions> &options, const std::shared_ptr<FatropPrinter> &printer);
        int eval_lag_hess(
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &lam) override;
        int eval_constr_jac(
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars) override;
        int solve_pd_sys(
            const double inertia_correction_w,
            const double inertia_correction_c,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_s,
            const FatropVecBF &sigma_total,
            const FatropVecBF &gradb_total) override;
        int solve_soc_rhs(
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_s,
            const FatropVecBF &constraint_violation) override;
        int compute_scalings(
            double &obj_scale,
            FatropVecBF &x_scales,
            FatropVecBF &lam_scales,
            const FatropVecBF &grad_curr);
        int eval_constraint_viol(
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            FatropVecBF &constraint_violation) override;
        int eval_obj_grad(
            double obj_scale,
            const FatropVecBF &primal_vars,
            FatropVecBF &gradient) override;
        int eval_obj(
            double obj_scale,
            const FatropVecBF &primal_vars,
            double &res) override;
        int eval_dual_inf(
            double obj_scale,
            const FatropVecBF &lam,
            const FatropVecBF &grad,
            FatropVecBF &du_inf) override;
        int initialize_slacks(
            FatropVecBF &s_curr) override;
        int initialize_dual(
            const FatropVecBF &grad,
            FatropVecBF &dlam,
            const FatropVecBF &zL,
            const FatropVecBF &zU) override;

        int get_bounds(
            FatropVecBF &lower,
            FatropVecBF &upper) const override
        {
            return ocp_->get_bounds(lower, upper);
        };
        int get_initial_sol_guess(
            FatropVecBF &initial) const override
        {
            return ocp_->get_initial_sol_guess(initial);
        };
        // int GetDefaultParams(
        //     FatropOptions &params) const override
        //     {
        //        return ocp_->GetDefaultParams(params);
        //     };
        NLPDims get_nlp_dims() const override;
        void finalize() override;
        void reset() override;

    public:
        std::shared_ptr<OCP> ocp_;
        OCPDims dims_;
        NLPDims nlpdims_;
        std::shared_ptr<OCPLinearSolver> ls_;
        std::shared_ptr<OCPScalingMethod> scaler_;
        std::shared_ptr<FatropOptions> options_;
        std::shared_ptr<FatropPrinter> printer_;
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
        bool it_ref;
    };
} // namespace fatrop
#endif //  OCPALGINCLUDED