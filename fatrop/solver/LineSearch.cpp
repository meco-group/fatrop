/*
 * Fatrop - A fast trajectory optimization solver
 * Copyright (C) 2022, 2023 Lander Vanroye, KU Leuven. All rights reserved.
 *
 * This file is part of Fatrop.
 *
 * Fatrop is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fatrop is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Fatrop.  If not, see <http://www.gnu.org/licenses/>. */
#include "solver/LineSearch.hpp"
using namespace fatrop;
using namespace std;
LineSearch::LineSearch(
    const shared_ptr<FatropOptions> &fatropparams,
    const shared_ptr<FatropNLP> &nlp,
    const shared_ptr<FatropData> &fatropdata, const std::shared_ptr<FatropPrinter> &printer) : AlgStrategy(fatropparams),
                                                                                               fatropnlp_(nlp),
                                                                                               fatropdata_(fatropdata), printer_(printer){};
inline fatrop_int LineSearch::eval_constr_viol_trial()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    fatrop_int res = fatropnlp_->eval_constraint_viol(
        fatropdata_->x_next,
        fatropdata_->s_next,
        fatropdata_->g_next);
    eval_cv_time += blasfeo_toc(&timer);
    eval_cv_count++;
    return res;
}
double LineSearch::eval_obj_trial()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    double res = 0.0;
    fatropnlp_->eval_obj(
        fatropdata_->obj_scale,
        fatropdata_->x_next,
        res);
    eval_obj_time += blasfeo_toc(&timer);
    eval_obj_count++;
    return res;
}
void LineSearch::reset()
{
    eval_cv_count = 0;
    eval_obj_count = 0;
    eval_cv_time = 0.;
    eval_obj_time = 0.;
}
fatrop_int LineSearch::update_trial_step(double alpha_pr, double alpha_du) const
{
    return fatropdata_->update_trial_step(alpha_pr, alpha_du);
};
fatrop_int LineSearch::initialize_second_order_correction() const
{
    // backup delta_x, delta_s and lam_calc
    fatropdata_->lam_calc_backup_ls.copy(fatropdata_->lam_calc);
    fatropdata_->delta_x_backup_ls.copy(fatropdata_->delta_x);
    fatropdata_->delta_s_backup_ls.copy(fatropdata_->delta_s);
    fatropdata_->g_soc.copy(fatropdata_->g_curr);
    return 0;
};
fatrop_int LineSearch::exit_second_order_correction() const
{
    // restore delta_x, delta_s and lam_calc
    fatropdata_->lam_calc.copy(fatropdata_->lam_calc_backup_ls);
    fatropdata_->delta_x.copy(fatropdata_->delta_x_backup_ls);
    fatropdata_->delta_s.copy(fatropdata_->delta_s_backup_ls);
    return 0;
};
fatrop_int LineSearch::compute_second_order_correction(double alpha) const
{
    axpy(alpha, fatropdata_->g_soc, fatropdata_->g_next, fatropdata_->g_soc);
    fatrop_int res = fatropnlp_->solve_soc_rhs(fatropdata_->delta_x, fatropdata_->lam_calc, fatropdata_->delta_s, fatropdata_->g_soc);
    if (res != 0)
    {
        printer_->level(1) << "SolveSOC failed" << endl;
    }
    return res;
};

BackTrackingLineSearch::BackTrackingLineSearch(
    const shared_ptr<FatropOptions> &fatropparams,
    const shared_ptr<FatropNLP> &nlp,
    const shared_ptr<FatropData> &fatropdata,
    const shared_ptr<Filter> &filter,
    const shared_ptr<Journaller> &journaller,
    const shared_ptr<FatropPrinter> &printer)
    : LineSearch(fatropparams, nlp, fatropdata, printer), filter_(filter), journaller_(journaller)
{
    initialize();
    fatrop_params_->register_option(BooleanOption("accept_every_trial_step", "accept every trial step", &accept_every_trial_step, false));
    fatrop_params_->register_option(NumericOption::lower_bounded("s_phi", "s_phi", &s_phi, 2.3, 0.0));
    fatrop_params_->register_option(NumericOption::lower_bounded("delta", "delta", &delta, 1.0, 0.0));
    fatrop_params_->register_option(NumericOption::lower_bounded("s_theta", "s_theta", &s_theta, 1.1, 0.0));
    fatrop_params_->register_option(NumericOption::lower_bounded("gamma_theta", "gamma_theta", &gamma_theta, 1e-12, 0.0));
    fatrop_params_->register_option(NumericOption::lower_bounded("gamma_phi", "gamma_phi", &gamma_phi, 1e-8, 0.0));
    fatrop_params_->register_option(NumericOption::lower_bounded("eta_phi", "eta_phi", &eta_phi, 1e-8, 0.0));
    fatrop_params_->register_option(NumericOption::lower_bounded("gamma_alpha", "gamma_alpha", &gamma_alpha, 0.05, 0.0));
    fatrop_params_->register_option(IntegerOption::lower_bounded("max_soc", "max_soc", &max_soc, 2, 0));
};
void BackTrackingLineSearch::initialize()
{
}
LineSearchInfo BackTrackingLineSearch::find_acceptable_trial_point(double mu, bool small_sd, bool from_backup)
{
    LineSearchInfo res;
    double alpha_primal = 1.0;
    double alpha_dual = 1.0;
    double alpha_primal_backup = 1.0;
    double alpha_dual_backup = 1.0;
    fatropdata_->maximum_step_size(alpha_primal, alpha_dual, MAX(1 - mu, 0.99));
    update_trial_step(alpha_primal, alpha_dual);
    double cv_curr = from_backup ? fatropdata_->constr_viol_sum_backup() : fatropdata_->constr_viol_sum_curr();
    double obj_curr = from_backup ? fatropdata_->obj_backup : fatropdata_->obj_curr;
    double barrier_curr = from_backup ? fatropdata_->eval_barrier_func_backup(mu) : fatropdata_->eval_barrier_func_curr(mu);
    obj_curr += barrier_curr;
    double lin_decr_curr = from_backup ? fatropdata_->fo_decr_obj_backup() : fatropdata_->fo_decr_obj_curr();
    double barrier_decr_curr = from_backup ? fatropdata_->eval_barrier_fo_decr_backup(mu) : fatropdata_->eval_barrier_fo_decr_curr(mu);
    lin_decr_curr += barrier_decr_curr;
    double theta_min = fatropdata_->theta_min;
    // calculation of alpha_min
    double alpha_min = gamma_alpha *
                       (lin_decr_curr > 0 ? gamma_theta
                                          : MIN(gamma_theta,
                                                cv_curr < theta_min ? MIN(-gamma_phi * cv_curr / lin_decr_curr, delta * pow(cv_curr, s_theta) / (pow(-lin_decr_curr, s_phi)))
                                                                    : -gamma_phi * cv_curr / lin_decr_curr));

    // cout << "alpha_min " << alpha_min << endl;
    // cout << "cv " << cv_curr << endl;
    // cout << "obj " << obj_curr << endl;
    // cout << "lindecr " << lin_decr_curr << endl;
    // const fatrop_int max_soc = 2;
    bool soc_step = false;
    double cv_soc_old = cv_curr;
    fatrop_int p = 0;
    fatrop_int no_alpha_trials = 1;
    for (fatrop_int ll = 1; ll < 500; ll++)
    {
        update_trial_step(alpha_primal, alpha_dual);
        if (alpha_primal < alpha_min)
        {
            res.ls = 0;
            return res;
        }
        eval_constr_viol_trial();
        double cv_next = fatropdata_->constr_viol_sum_next();
        double obj_next = eval_obj_trial();
        double barrier_next = fatropdata_->eval_barrier_func_trial(mu);
        obj_next += barrier_next;

        // cout << "cv_next: " << cv_next;
        // cout << "  obj_next: " << obj_next << endl;
        double alpha_primal_accent = (soc_step ? alpha_primal_backup : alpha_primal);
        if (filter_->is_acceptable(FilterData(0, obj_next, cv_next)))
        {
            // cout << filter_->GetSize() << endl;
            bool switch_cond = (lin_decr_curr < 0) && (alpha_primal_accent * pow(-lin_decr_curr, s_phi) > delta * pow(cv_curr, s_theta));
            bool armijo = obj_next - obj_curr < eta_phi * alpha_primal_accent * lin_decr_curr;
            if (switch_cond && (cv_curr <= theta_min))
            {
                // f-step
                if (armijo)
                {
                    (journaller_->it_curr).type = soc_step ? 'F' : 'f';
                    fatropdata_->accept_trial_step();
                    journaller_->it_curr.alpha_pr = alpha_primal;
                    journaller_->it_curr.alpha_du = alpha_dual;
                    res.ls = no_alpha_trials;
                    return res;
                }
            }
            else
            {
                // h-step
                // check sufficient decrease wrt current iterate
                if ((cv_next < (1.0 - gamma_theta) * cv_curr) || (obj_next < obj_curr - gamma_phi * cv_curr))
                {
                    if (!switch_cond || !(armijo))
                    {
                        filter_->augment(FilterData(0, obj_curr - gamma_phi * cv_curr, cv_curr * (1 - gamma_theta)));
                    }
                    (journaller_->it_curr).type = soc_step ? 'H' : 'h';
                    fatropdata_->accept_trial_step();
                    journaller_->it_curr.alpha_pr = alpha_primal;
                    journaller_->it_curr.alpha_du = alpha_dual;
                    res.ls = no_alpha_trials;
                    return res;
                }
            }
            res.last_rejected_by_filter = false;
        }
        else
        {
            res.last_rejected_by_filter = true;
            if (soc_step)
            {
                // abort soc
                p = max_soc;
            }
            if (ll == 1)
            {
                res.first_rejected_by_filter = true;
            }
        }
        // todo change iteration number from zero to real iteration number
        if (small_sd)
        {
            (journaller_->it_curr).type = 's';
            fatropdata_->accept_trial_step();
            journaller_->it_curr.alpha_pr = alpha_primal;
            journaller_->it_curr.alpha_du = alpha_dual;
            res.ls = -1;
            return res;
        }
        if (accept_every_trial_step)
        {
            (journaller_->it_curr).type = 'a';
            fatropdata_->accept_trial_step();
            journaller_->it_curr.alpha_pr = alpha_primal;
            journaller_->it_curr.alpha_du = alpha_dual;
            res.ls = 1;
            return res;
        }
        if (soc_step && (p >= max_soc || (ll > 1 && (cv_next > 0.99 * cv_soc_old) && p > 1)))
        {
            // deactivate soc
            soc_step = false;
            alpha_primal = alpha_primal_backup;
            alpha_dual = alpha_dual_backup;
            exit_second_order_correction();
            // todo cache these variables
            fatropdata_->compute_delta_z();
        }
        if (!soc_step && (ll == 1 && max_soc > 0) && cv_next > cv_curr)
        {
            // activate soc
            // cout << "trying soc " << endl;
            soc_step = true;
            alpha_primal_backup = alpha_primal;
            alpha_dual_backup = alpha_dual;
            initialize_second_order_correction();
        }
        if (soc_step)
        {
            if (ll > 1)
            {
                cv_soc_old = cv_next;
            }
            fatrop_int res = compute_second_order_correction(alpha_primal);
            if (res == 0)
            {
                // cout << "size of soc x step " << L1(fatropdata_->delta_x) << endl;
                fatropdata_->compute_delta_z();
                fatropdata_->maximum_step_size(alpha_primal, alpha_dual, MAX(1 - mu, 0.99));
                // cout << "alpha_primal " << alpha_primal << endl;
                // cout << "alpha_dual " << alpha_dual << endl;
                // cout << "soc" << endl;
                p = p + 1;
            }
            else
            {
                soc_step = false;
                alpha_primal = alpha_primal_backup;
                alpha_dual = alpha_dual_backup;
                exit_second_order_correction();
                // todo cache these variables
                fatropdata_->compute_delta_z();
            }
        }
        else
        {
            // cut back alpha_primal
            no_alpha_trials++;
            alpha_primal *= 0.50;
        }
    }
    assert(false);
    res.ls = 0;
    return res;
};