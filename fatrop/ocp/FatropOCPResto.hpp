#pragma once
#include "fatrop/templates/NLPAlg.hpp"
namespace fatrop
{
    class FatropOCPResto : public FatropNLP
    {
        FatropOCPResto(const std::shared_ptr<FatropOCP> &orig) : orig_(orig), orig_dims_(orig->get_nlp_dims()), lower_(orig_dims_.nineqs), upper_(orig_dims_.nineqs), upper_bounded_(orig_dims_.nineqs), lower_bounded_(orig_dims_.nineqs), slack_dummy_(orig_dims_.nineqs)
        {
            auto lower_v = lower_[0];
            auto upper_v = upper_[0];
            orig_->get_bounds(lower_v, upper_v);
            n_p = orig_dims_.nineqs;
            n_n = orig_dims_.nineqs;
            for (fatrop_int i = 0; i < orig_dims_.nineqs; ++i)
            {
                bool lower_bounded = !std::isinf(lower_v.at(i));
                bool upper_bounded = !std::isinf(upper_v.at(i));
                upper_bounded_[i] = upper_bounded;
                lower_bounded_[i] = lower_bounded;
            }
        }
        virtual fatrop_int eval_lag_hess(
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            const FatropVecBF &lam) override
        {
            orig_->eval_lag_hess(obj_scale, primal_vars, slack_vars, lam);
            return 0;
        } ;

        virtual fatrop_int eval_constr_jac(
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars) override
        {
            FatropVecBF slack_dummy_v = slack_dummy_[0];
            update_slack_vars(slack_vars, slack_dummy_v);
            orig_->eval_constr_jac(primal_vars, slack_dummy_v);
            return 0;
        };
        void update_slack_vars(const FatropVecBF &slack_vars, FatropVecBF &slack_dummy)
        {
            fatrop_int offs_n = orig_dims_.nineqs;
            fatrop_int offs_p = orig_dims_.nineqs + n_n;
            for (fatrop_int i = 0; i < orig_dims_.nineqs; ++i)
            {
                slack_dummy.at(i) = slack_vars.at(i) - slack_vars.at(i + offs_n) + slack_vars.at(i + offs_p);
            }
        }
        virtual fatrop_int eval_constraint_viol(
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            FatropVecBF &constraint_violation) override
        {
            FatropVecBF slack_dummy_v = slack_dummy_[0];
            update_slack_vars(slack_vars, slack_dummy_v);
            orig_->eval_constraint_viol(primal_vars, slack_dummy_v, constraint_violation);
            return 0;
        };
        virtual fatrop_int eval_obj_grad(
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            FatropVecBF &gradient_x,
            FatropVecBF &gradient_s)
        {
            orig_->eval_obj_grad(obj_scale, primal_vars, slack_vars, gradient_x, gradient_s);
            gradient_x.block(orig_dims_.nvars, n_n + n_p) = rho;
            return 0;
        };
        virtual fatrop_int eval_obj(
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            double &res) override
        {
            res = orig_->eval_obj(obj_scale, primal_vars, slack_vars, res);
            res += rho * sum(primal_vars.block(orig_dims_.nvars, n_n + n_p));
            return 0;
        };
        void get_initialization(const FatropVecBF &primal_vars_orig, FatropVecBF &intialization)
        {
            intialization.block(0, orig_dims_.nvars).copy(primal_vars_orig.block(0, orig_dims_.nvars));
        }
        virtual fatrop_int eval_ineqs(
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            FatropVecBF &ineqs){
                return 0;
        };
        virtual fatrop_int eval_dual_inf(
            double obj_scale,
            const FatropVecBF &lam,
            const FatropVecBF &grad_x,
            const FatropVecBF &grad_s,
            FatropVecBF &du_inf) = 0;
        virtual fatrop_int solve_pd_sys(
            const double inertia_correction_w,
            const double inertia_correction_c,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_s,
            const FatropVecBF &sigma_total,
            const FatropVecBF &gradb_total) = 0;
        virtual fatrop_int solve_soc_rhs(
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_s,
            const FatropVecBF &cosntraint_violation) = 0;
        virtual NLPDims get_nlp_dims() const = 0;
        virtual fatrop_int compute_scalings(
            double &obj_scale,
            FatropVecBF &x_scales,
            FatropVecBF &lam_scales,
            const FatropVecBF &grad_curr) = 0;
        virtual fatrop_int initialize_slacks(
            FatropVecBF &s_curr) = 0;
        virtual fatrop_int initialize_dual(
            const FatropVecBF &grad,
            FatropVecBF &dlam,
            const FatropVecBF &zL,
            const FatropVecBF &zU) = 0;
        virtual fatrop_int get_bounds(
            FatropVecBF &lower,
            FatropVecBF &upper) const = 0;
        virtual fatrop_int get_initial_sol_guess(
            FatropVecBF &initial) const = 0;
        std::shared_ptr<FatropOCP> orig_;
        NLPDims orig_dims_;
        NLPDims this_dims_;
        FatropMemoryVecBF lower_, upper_;
        std::vector<bool> upper_bounded_;
        std::vector<bool> lower_bounded_;
        fatrop_int n_p = 0;
        fatrop_int n_n = 0;
        FatropMemoryVecBF slack_dummy_;
        double rho = 1e4;
    };
};