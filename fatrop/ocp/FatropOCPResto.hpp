#pragma once
#include "fatrop/templates/NLPAlg.hpp"
namespace fatrop
{
    class FatropOCPResto : public FatropNLP
    {
        FatropOCPResto(const std::shared_ptr<FatropOCP> &orig) : orig_(orig), orig_dims_(orig->get_nlp_dims()), lower_(orig_dims_.nineqs), upper_(orig_dims_.nineqs), upper_bounded_(orig_dims_.nineqs), lower_bounded_(orig_dims_.nineqs), slack_dummy_(orig_dims_.nineqs), sigma_dummy_(orig_dims_.nineqs), gradb_dummy_(orig_dims_.nineqs)
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
            this_dims_.nvars = orig_dims_.nvars;
            this_dims_.nineqs = orig_dims_.nineqs + n_n + n_p;
            this_dims_.neqs = orig_dims_.neqs;
        }
        virtual fatrop_int eval_lag_hess(
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            const FatropVecBF &lam) override
        {
            orig_->eval_lag_hess(obj_scale, primal_vars, slack_vars, lam);
            return 0;
        };

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
            FatropVecBF &gradient_s) override
        {
            orig_->eval_obj_grad(obj_scale, primal_vars, slack_vars, gradient_x, gradient_s);
            gradient_s.block(orig_dims_.nineqs, n_n + n_p) = rho;
            return 0;
        };
        virtual fatrop_int eval_obj(
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            double &res) override
        {
            res = orig_->eval_obj(obj_scale, primal_vars, slack_vars, res);
            res += rho * sum(slack_vars.block(orig_dims_.nineqs, n_n + n_p));
            return 0;
        };
        void get_initialization(const FatropVecBF &primal_vars_orig, FatropVecBF &intialization)
        {
            intialization.block(0, orig_dims_.nvars).copy(primal_vars_orig.block(0, orig_dims_.nvars));
        }
        virtual fatrop_int eval_dual_inf(
            double obj_scale,
            const FatropVecBF &lam,
            const FatropVecBF &grad_x,
            const FatropVecBF &grad_s,
            FatropVecBF &du_inf_x, FatropVecBF &du_inf_s_wo_z) override
        {
            orig_->eval_dual_inf(obj_scale, lam, grad_x, grad_s, du_inf_x, du_inf_s_wo_z);
            auto lam_I = lam.block(orig_dims_.neqs - orig_dims_.nineqs, orig_dims_.nineqs);
            axpby(-1.0, lam_I, 0.0, lam_I, du_inf_s_wo_z.block(orig_dims_.nineqs, n_n));
            axpby(1.0, lam_I, 0.0, lam_I, du_inf_s_wo_z.block(orig_dims_.nineqs + n_n, n_p));
        }
        void update_sigma_gradb(double inertia, const FatropVecBF &sigma_s, const FatropVecBF &sigma_n, const FatropVecBF &sigma_p, const FatropVecBF &gradb_s, const FatropVecBF &gradb_n, const FatropVecBF &gradb_p, const FatropVecBF &sigma_update, const FatropVecBF &gradb_update)
        {
            for (int i = 0; i < orig_dims_.nineqs; ++i)
            {
                double sigma_updt = 1.0 / (1.0 / (sigma_s.at(i) + inertia) + 1.0 / (sigma_n.at(i) + inertia) + 1.0 / (sigma_p.at(i) + inertia));
                sigma_update.at(i) = sigma_updt - inertia;
                gradb_update.at(i) = (gradb_s.at(i) / (sigma_s.at(i) + inertia) - gradb_n.at(i) / (sigma_n.at(i) + inertia) + gradb_p.at(i) / (sigma_p.at(i) + inertia)) * sigma_updt;
            }
        }
        void update_delta_snp(double inertia, const FatropVecBF &sigma_s, const FatropVecBF &sigma_n, const FatropVecBF &sigma_p, const FatropVecBF &gradb_s, const FatropVecBF &gradb_n, const FatropVecBF &gradb_p, const FatropVecBF &lam_I, const FatropVecBF &delta_s, const FatropVecBF &delta_n, const FatropVecBF &delta_p)
        {
            for (int i = 0; i < orig_dims_.nineqs; ++i)
            {
                double lam_I_i = lam_I.at(i);
                delta_s.at(i) = (-gradb_s.at(i) + lam_I_i) / (sigma_s.at(i) + inertia);
                delta_n.at(i) = (-gradb_n.at(i) - lam_I_i) / (sigma_n.at(i) + inertia);
                delta_p.at(i) = (-gradb_p.at(i) + lam_I_i) / (sigma_p.at(i) + inertia);
            }
        }

        virtual fatrop_int solve_pd_sys(
            const double inertia_correction_w,
            const double inertia_correction_c,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_s,
            const FatropVecBF &sigma_total,
            const FatropVecBF &gradb_total) override
        {
            inertia_correction_w_cache = inertia_correction_w;
            auto lam_I = lam.block(orig_dims_.neqs - orig_dims_.nineqs, orig_dims_.nineqs);
            auto sigma_s = sigma_total.block(0, orig_dims_.nineqs);
            auto sigma_n = sigma_total.block(orig_dims_.nineqs, n_n);
            auto sigma_p = sigma_total.block(orig_dims_.nineqs + n_n, n_p);
            auto gradb_s = gradb_total.block(0, orig_dims_.nineqs);
            auto gradb_n = gradb_total.block(orig_dims_.nineqs, n_n);
            auto gradb_p = gradb_total.block(orig_dims_.nineqs + n_n, n_p);
            auto delta_s_or = delta_s.block(0, orig_dims_.nineqs);
            auto delta_n = delta_s.block(orig_dims_.nineqs, n_n);
            auto delta_p = delta_s.block(orig_dims_.nineqs + n_n, n_p);
            update_sigma_gradb(inertia_correction_w, sigma_s, sigma_n, sigma_p, gradb_s, gradb_n, gradb_p, sigma_dummy_[0], gradb_dummy_[0]);
            int ret = orig_->solve_pd_sys(inertia_correction_w, inertia_correction_c, ux, lam, delta_s, sigma_dummy_[0], gradb_dummy_[0]);
            update_delta_snp(inertia_correction_w, sigma_s, sigma_n, sigma_p, gradb_s, gradb_n, gradb_p, lam_I, delta_s_or, delta_n, delta_p);
            return ret;
        };
        virtual fatrop_int solve_soc_rhs(
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_s,
            const FatropVecBF &cosntraint_violation) override
        {
            auto lam_I = lam.block(orig_dims_.neqs - orig_dims_.nineqs, orig_dims_.nineqs);
            auto sigma_s = sigma_dummy_[0].block(0, orig_dims_.nineqs);
            auto sigma_n = sigma_dummy_[0].block(orig_dims_.nineqs, n_n);
            auto sigma_p = sigma_dummy_[0].block(orig_dims_.nineqs + n_n, n_p);
            auto gradb_s = gradb_dummy_[0].block(0, orig_dims_.nineqs);
            auto gradb_n = gradb_dummy_[0].block(orig_dims_.nineqs, n_n);
            auto gradb_p = gradb_dummy_[0].block(orig_dims_.nineqs + n_n, n_p);
            auto delta_s_or = delta_s.block(0, orig_dims_.nineqs);
            auto delta_n = delta_s.block(orig_dims_.nineqs, n_n);
            auto delta_p = delta_s.block(orig_dims_.nineqs + n_n, n_p);
            int ret = orig_->solve_soc_rhs(ux, lam, delta_s, cosntraint_violation);
            update_delta_snp(inertia_correction_w_cache, sigma_s, sigma_n, sigma_p, gradb_s, gradb_n, gradb_p, lam_I, delta_s_or, delta_n, delta_p);
            return ret;
        }
        virtual NLPDims get_nlp_dims() const override { return this_dims_; };
        virtual fatrop_int compute_scalings(
            double &obj_scale,
            FatropVecBF &x_scales,
            FatropVecBF &lam_scales,
            const FatropVecBF &grad_curr_x, const FatropVecBF &grad_curr_s) override
        {
            return orig_->compute_scalings(obj_scale, x_scales, lam_scales, grad_curr_x, grad_curr_s);
        };
        virtual fatrop_int initialize_slacks(
            FatropVecBF &s_curr) override
        {
            auto s_curr_or = s_curr.block(0, orig_dims_.nineqs);
            auto n_curr = s_curr.block(orig_dims_.nineqs, n_n);
            auto p_curr = s_curr.block(orig_dims_.nineqs + n_n, n_p);
            orig_->initialize_slacks(s_curr_or);
            // set zero for now TODO use quadratic formula here
            n_curr = 0.0;
            p_curr = 0.0;
        }
        virtual fatrop_int initialize_dual(
            const FatropVecBF &grad_x,
            const FatropVecBF &grad_s,
            FatropVecBF &dlam,
            const FatropVecBF &zL,
            const FatropVecBF &zU) override = 0;
        virtual fatrop_int get_bounds(
            FatropVecBF &lower,
            FatropVecBF &upper) const override = 0;
        virtual fatrop_int get_initial_sol_guess(
            FatropVecBF &initial) const override = 0;
        std::shared_ptr<FatropOCP> orig_;
        NLPDims orig_dims_;
        NLPDims this_dims_;
        FatropMemoryVecBF lower_, upper_;
        std::vector<bool> upper_bounded_;
        std::vector<bool> lower_bounded_;
        fatrop_int n_p = 0;
        fatrop_int n_n = 0;
        FatropMemoryVecBF slack_dummy_;
        FatropMemoryVecBF sigma_dummy_;
        FatropMemoryVecBF gradb_dummy_;
        double inertia_correction_w_cache = 0.0;
        double rho = 1e4;
    };
};