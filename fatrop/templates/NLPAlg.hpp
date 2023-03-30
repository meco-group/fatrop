#ifndef NLPINCLUDED
#define NLPINCLUDED
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
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
        virtual int eval_lag_hess(
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &lam) = 0;
        virtual int eval_constr_jac(
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars) = 0;
        virtual int eval_constraint_viol(
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            FatropVecBF &constraint_violation) = 0;
        virtual int eval_obj_grad(
            double obj_scale,
            const FatropVecBF &primal_vars,
            FatropVecBF &gradient) = 0;
        virtual int eval_obj(
            double obj_scale,
            const FatropVecBF &primal_vars,
            double &res) = 0;
        virtual int eval_dual_inf(
            double obj_scale,
            const FatropVecBF &lam,
            const FatropVecBF &grad,
            FatropVecBF &du_inf) = 0;
        virtual int solve_pd_sys(
            const double inertia_correction_w,
            const double inertia_correction_c,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_s,
            const FatropVecBF &sigma_total,
            const FatropVecBF &gradb_total) = 0;
        virtual int solve_soc_rhs(
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_s,
            const FatropVecBF &cosntraint_violation) = 0;
        virtual NLPDims get_nlp_dims() const = 0;
        virtual int compute_scalings(
            double &obj_scale,
            FatropVecBF &x_scales,
            FatropVecBF &lam_scales,
            const FatropVecBF &grad_curr) = 0;
        virtual int initialize_slacks(
            FatropVecBF &s_curr) = 0;
        virtual int initialize_dual(
            const FatropVecBF &grad,
            FatropVecBF &dlam,
            const FatropVecBF &zL,
            const FatropVecBF &zU) = 0;
        virtual int get_bounds(
            FatropVecBF &lower,
            FatropVecBF &upper) const = 0;
        virtual int get_initial_sol_guess(
            FatropVecBF &initial) const = 0;
        // virtual int GetDefaultParams(
        //     FatropOptions &params) const = 0;
        virtual int Callback(FatropVecBF& primal_vars){return 0;};
        virtual void finalize(){};
        virtual void reset(){};
    };
} // namespace fatrop
#endif // NLPINCLUDED