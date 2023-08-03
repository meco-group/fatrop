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
#include "solver/FatropAlg.hpp"
using namespace fatrop;
using namespace std;

FatropAlg::FatropAlg(
    const shared_ptr<FatropNLP> &fatropnlp,
    const shared_ptr<FatropData> &fatropdata,
    const shared_ptr<FatropOptions> &fatropparams,
    const shared_ptr<Filter> &filter,
    const shared_ptr<LineSearch> &linesearch,
    const shared_ptr<Journaller> &journaller,
    const shared_ptr<FatropPrinter> &printer)
    : fatropnlp_(fatropnlp),
      fatropdata_(fatropdata),
      fatropoptions_(fatropparams),
      filter_(filter),
      linesearch_(linesearch),
      journaller_(journaller),
      printer_(printer)
{

    fatropoptions_->register_option(NumericOption::lower_bounded("tol", "tolerance", &tol, 1e-8, 0.0));
    fatropoptions_->register_option(NumericOption::lower_bounded("acceptable_tol", "acceptable tolerance", &acceptable_tol, 1e-6, 0.0));
    fatropoptions_->register_option(IntegerOption::lower_bounded("max_watchdog_steps", "maximum number of watchdog steps", &max_watchdog_steps, 4, 0));
    fatropoptions_->register_option(IntegerOption::lower_bounded("acceptable_iter", "acceptable iter", &acceptable_iter, 15, 0));
    fatropoptions_->register_option(NumericOption::lower_bounded("lammax", "lammax", &lammax, 1e3, 0.0));
    fatropoptions_->register_option(NumericOption::lower_bounded("mu_init", "mu_init", &mu0, 1e2, 0.0));
    fatropoptions_->register_option(NumericOption::lower_bounded("kappa_eta", "kappa_eta", &kappa_eta, 10.0, 0.0));
    fatropoptions_->register_option(NumericOption::lower_bounded("kappa_mu", "kappa_mu", &kappa_mu, 0.2, 0.0));
    fatropoptions_->register_option(NumericOption::lower_bounded("theta_mu", "theta_mu", &theta_mu, 1.5, 0.0));
    fatropoptions_->register_option(NumericOption::lower_bounded("delta_w0", "delta_w0", &delta_w0, 1e-4, 0.0));
    fatropoptions_->register_option(NumericOption::lower_bounded("delta_wmin", "delta_wmin", &delta_wmin, 1e-20, 0.0));
    fatropoptions_->register_option(NumericOption::lower_bounded("kappa_wmin", "kappa_wmin", &kappa_wmin, 1.0 / 3.0, 0.0));
    fatropoptions_->register_option(NumericOption::lower_bounded("kappa_wplus", "kappa_wplus", &kappa_wplus, 8.0, 0.0));
    fatropoptions_->register_option(NumericOption::lower_bounded("kappa_wplusem", "kappa_wplusem", &kappa_wplusem, 100.0, 0.0));
    fatropoptions_->register_option(NumericOption::lower_bounded("delta_c_stripe", "delta_c_stripe", &delta_c_stripe, 1e-2, 0.0));
    fatropoptions_->register_option(NumericOption::lower_bounded("kappa_c", "kappa_c", &kappa_c, 0.25, 0.0));
    fatropoptions_->register_option(BooleanOption("warm_start_init_point", "warm_start_init_point", &warm_start_init_point, false));
    fatropoptions_->register_option(NumericOption::lower_bounded("theta_min", "theta_min", &theta_min, 1e-4, 0.0));
    fatropoptions_->register_option(BooleanOption("recalc_y", "recalc_y", &recalc_y, false));
    fatropoptions_->register_option(NumericOption::lower_bounded("recalc_y_feas_tol", "recalc_y_feas_tol", &recalc_y_feas_tol, 1e-6, 0.0));
    initialize();
    fatropnlp_->get_initial_sol_guess(fatropdata_->x_initial);
    fatropnlp->get_bounds(fatropdata->s_lower_orig, fatropdata->s_upper_orig);
    fatropdata->relax_bounds();
}
void FatropAlg::initialize()
{
    maxiter = fatropoptions_->maxiter;
    kappa_d = fatropoptions_->kappa_d;
    fatropdata_->initialize();
    linesearch_->initialize();
}
void FatropAlg::reset()
{
    filter_->reset();
    fatropdata_->reset();
    journaller_->reset();
    fatropnlp_->reset();
    linesearch_->reset();
    stats = FatropStats();
}
void FatropAlg::set_bounds(const vector<double> &lower, const vector<double> &upper)
{
    fatropdata_->s_lower_orig = lower;
    fatropdata_->s_upper_orig = upper;
};
void FatropAlg::set_initial(const vector<double> &initial)
{
    fatropdata_->x_initial = initial;
};
void FatropAlg::get_solution(vector<double> &sol)
{
    fatropdata_->x_curr.copyto(sol);
};
fatrop_int FatropAlg::optimize()
{
    // bool first_try_watchdog = this->first_try_watchdog;
    fatrop_int no_watch_dog_steps_taken = 0;
    fatrop_int max_watchdog_steps = this->max_watchdog_steps;
    initialize();
    fatrop_int no_conse_small_sd = false;
    fatrop_int filter_reseted = 0;
    fatrop_int no_no_full_steps = 0;
    fatrop_int no_no_full_steps_bc_filter = 0;
    fatrop_int no_acceptable_steps = 0;
    // double delta_w_last_backup = 0.;
    bool restore_watchdog_step = false;
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    reset();
    const double mu_min = tol / 10;
    double mu = mu0;
    double delta_w_last = 0.0;
    LineSearchInfo lsinfo;
    eval_constr_jac(); // todo twice evaluation
    eval_obj_grad_curr();
    if (warm_start_init_point)
    {
        fatropnlp_->initialize_slacks(
            fatropdata_->s_curr);
        fatropdata_->warmstart_dual();
        fatropdata_->bound_z();
    }
    else
    {
        fatrop_int initialization_res = perform_initializiation();
        if (initialization_res == 0 && fatropdata_->delta_dual_max() < lammax)
        {
            printer_->level(1) << "accepted lam " << endl;
            fatropdata_->accept_dual_initializiaton();
        }
        else
        {
            printer_->level(1) << "rejected lam " << endl;
            fatropdata_->lam_curr.SetConstant(0.0);
        }
    }
    fatropdata_->bound_slacks();
    eval_constr_viol_curr();
    fatropdata_->theta_min = theta_min * MAX(1.0, fatropdata_->constr_viol_sum_curr());
    double theta_max = 1e4 * fatropdata_->constr_viol_sum_curr();
    filter_->augment(FilterData(0, std::numeric_limits<double>::infinity(), theta_max));
    fatrop_int ls = 0;
    double deltaw = 0;
    double deltac = 0.0;
    bool watch_dog_step = false;
    for (fatrop_int i = 0; i < maxiter; i++)
    {
        fatropdata_->obj_curr = eval_objective_curr();
        // if (fatropdata_->LamLinfCurr() > 1e12)
        // {
        //     cout << "huge Lagrange multipliers -> set to zero" << endl;
        //     fatropdata_->lam_curr.SetConstant(0.0);
        // }
        eval_constr_jac();    // needed for dual inf
        eval_obj_grad_curr(); // needed for dual inf
        eval_dual_infeasiblity();
        IterationData &it_curr = journaller_->it_curr;
        it_curr.iter = i;
        it_curr.mu = mu;
        it_curr.objective = eval_objective_curr();
        it_curr.constraint_violation = fatropdata_->constr_viol_max_curr();
        it_curr.du_inf = fatropdata_->dual_inf_max_curr();
        it_curr.ls = ls;
        it_curr.reg = deltaw;
        if (no_no_full_steps_bc_filter >= 5)
        {
            bool reset_filter = (filter_reseted <= 5);
            if (reset_filter)
            {
                printer_->level(1) << "resetted filter " << endl;
                filter_reseted++;
                filter_->reset();
                no_no_full_steps_bc_filter = 0;
                filter_->augment(FilterData(0, std::numeric_limits<double>::infinity(), theta_max));
            }
        }
        if (no_no_full_steps >= 10)
        {
            if (((max_watchdog_steps > 0) && !watch_dog_step))
            {
                // activate watchdog procedure
                // backup x_k
                fatropdata_->backup_curr();
                // delta_w_last_backup = delta_w_last;
                watch_dog_step = true;
                // no_no_full_steps = 0;
                no_watch_dog_steps_taken = 0;
            }
        }
        journaller_->push();
        journaller_->print_iterations();
        double emu = fatropdata_->e_mu_curr(0.0);
        if (emu < acceptable_tol)
        {
            no_acceptable_steps++;
        }
        else
        {
            no_acceptable_steps = 0;
        }

        if (emu < tol || (no_acceptable_steps >= acceptable_iter) || ((no_conse_small_sd == 2) && (mu <= mu_min)))
        // if (emu < tol)
        {
            double total_time = blasfeo_toc(&timer);
            journaller_->print_iterations();
            if (no_conse_small_sd == 2)
            {
                printer_->level(1) << "WARNING fatrop returned bc of very small search direction" << endl;
            }
            if (emu > tol && no_acceptable_steps >= acceptable_iter)
            {
                printer_->level(1) << "WARNING fatrop returned acceptable tolerance" << endl;
            }
            printer_->level(1) << "found solution :) " << endl;
            stats.eval_cv_count += linesearch_->eval_cv_count;
            stats.eval_obj_count += linesearch_->eval_obj_count;
            stats.eval_cv_time += linesearch_->eval_cv_time;
            stats.eval_obj_time += linesearch_->eval_obj_time;
            stats.time_total = total_time;
            stats.iterations_count = i;
            if (printer_->print_level() > 0)
            {
                stats.print();
            }
            fatropnlp_->finalize();
            return 0;
        }
        // update mu
        // todo make a seperate class
        while (!watch_dog_step && mu > mu_min && (fatropdata_->e_mu_curr(mu) <= kappa_eta * mu || (no_conse_small_sd == 2)))
        {
            mu = MAX(mu_min, MIN(kappa_mu * mu, pow(mu, theta_mu)));
            filter_reseted = 0;
            filter_->reset();
            no_no_full_steps_bc_filter = 0;
            filter_->augment(FilterData(0, std::numeric_limits<double>::infinity(), theta_max));
            if (no_conse_small_sd == 2)
            {
                // cout << "small search direction" << endl;
                no_conse_small_sd = 0;
                break;
            }
            no_no_full_steps = 0;
            // the following break statement prohibits 'fast' mu updates, at leat one iteration per mu update
            // break;
        }
        // Hessian is necessary for calculating search direction
        eval_lag_hess();
        // todo make an update class for regularization
        double deltac_candidate = delta_c_stripe * pow(mu, kappa_c);
        deltaw = 0.0;
        deltac = 0.0;
        fatropdata_->evaluate_barrier_quantities(mu);
        fatrop_int regularity = -1;
        fatrop_int increase_counter = 0;
        if (!restore_watchdog_step)
        {
            while (regularity != 0)
            {
                regularity = solve_pd_sys(deltaw, deltac, mu);
                if (deltac == 0 && regularity < 0)
                {
                    printer_->level(1) << "degenerate Jacobian" << endl;
                    deltac = deltac_candidate;
                }
                if (regularity > 0) // regularization is necessary
                {
                    if (increase_counter == 0)
                    {
                        deltaw = (delta_w_last == 0.0) ? delta_w0 : MAX(delta_wmin, kappa_wmin * delta_w_last);
                    }
                    else
                    {
                        deltaw = (delta_w_last == 0.0) ? kappa_wplusem * deltaw : kappa_wplus * deltaw;
                    }
                    increase_counter++;
                }
            }
            if (deltaw > 0.)
                delta_w_last = deltaw;
        }
        fatropdata_->compute_delta_z();
        // cout << "norm dzL " << Linf(fatropdata_->delta_zL) << endl;
        // cout << "norm dzU " << Linf(fatropdata_->delta_zU) << endl;
        // cout << "norm delta_s " << Linf(fatropdata_->delta_s) << endl;
        // cout << "norm delta_x " << Linf(fatropdata_->delta_x) << endl;
        // cout << "norm delta_lam " << Linf(fatropdata_->lam_calc) << endl;
        double stepsize = max(LinfScaled(fatropdata_->delta_x, fatropdata_->x_curr), LinfScaled(fatropdata_->delta_s, fatropdata_->s_curr));
        bool small_search_direction_curr = stepsize < 1e-14;
        lsinfo = linesearch_->find_acceptable_trial_point(mu, small_search_direction_curr || watch_dog_step, watch_dog_step);
        if (recalc_y && (deltac == 0.0) && (fatropdata_->constr_viol_max_curr() < recalc_y_feas_tol))
        {
            fatropnlp_->initialize_dual(
                fatropdata_->grad_curr,
                fatropdata_->lam_curr,
                fatropdata_->zL_curr,
                fatropdata_->zU_curr);
        }
        fatropdata_->relax_bounds_var(mu);
        fatropdata_->modify_dual_bounds(mu);
        ls = lsinfo.ls;
        if (watch_dog_step && no_watch_dog_steps_taken == 0)
        {
            fatropdata_->backup_delta();
        }
        if (watch_dog_step)
        {
            if (ls == 1)
            {
                // accept watchdog step -- continue
                printer_->level(1) << "accepted watchdog step" << endl;
                watch_dog_step = false;
            }
            else
            {
                no_watch_dog_steps_taken++;
                if (no_watch_dog_steps_taken >= max_watchdog_steps)
                {
                    // reject watchdog step -- go back to x_k
                    printer_->level(1) << "rejected watchdog step" << endl;
                    it_curr.type = 'x';
                    fatropdata_->restore_backup();
                    // delta_w_last = delta_w_last_backup;
                    watch_dog_step = false;
                    // todo make use of delta_x_backup and delta_s_backup
                    restore_watchdog_step = true;
                    i--;
                    continue;
                };
                it_curr.type = 'w';
            }
        }

        if (ls == 0)
        {
            printer_->level(1) << "error: restoration phase not implemented yet" << endl;
            return 1;
        }
        if (watch_dog_step || ls == 1)
        {
            no_no_full_steps = 0;
            no_no_full_steps_bc_filter = 0;
        }
        else
        {
            ++no_no_full_steps;
            if (lsinfo.last_rejected_by_filter)
            {
                ++no_no_full_steps_bc_filter;
            }
            else
            {
                no_no_full_steps_bc_filter = 0;
            }
        }
        if (small_search_direction_curr)
        {
            it_curr.type = 's';
            no_conse_small_sd++;
            if (!no_conse_small_sd)
            {
                // take the full step
                cout << "full step small sd " << endl;
            }
        }
        else
        {
            no_conse_small_sd = 0;
        }
        restore_watchdog_step = false;
        // if linesearch unsuccessful -> resto phase
    }
    journaller_->print_iterations();
    return 1;
}
fatrop_int FatropAlg::eval_lag_hess()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    fatrop_int res =
        fatropnlp_->eval_lag_hess(
            fatropdata_->obj_scale,
            fatropdata_->x_curr,
            fatropdata_->lam_curr);
    stats.eval_hess_time += blasfeo_toc(&timer);
    stats.eval_hess_count++;
    return res;
}
fatrop_int FatropAlg::eval_constr_jac()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    fatrop_int res = fatropnlp_->eval_constr_jac(
        fatropdata_->x_curr,
        fatropdata_->s_curr);
    stats.eval_jac_time += blasfeo_toc(&timer);
    stats.eval_jac_count++;
    return res;
}
fatrop_int FatropAlg::eval_constr_viol_curr()
{

    blasfeo_timer timer;
    blasfeo_tic(&timer);
    fatrop_int res = fatropnlp_->eval_constraint_viol(
        fatropdata_->x_curr,
        fatropdata_->s_curr,
        fatropdata_->g_curr);
    stats.eval_cv_time += blasfeo_toc(&timer);
    stats.eval_cv_count++;
    return res;
}
fatrop_int FatropAlg::eval_obj_grad_curr()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    fatrop_int res = fatropnlp_->eval_obj_grad(
        fatropdata_->obj_scale,
        fatropdata_->x_curr,
        fatropdata_->grad_curr);
    stats.eval_grad_time += blasfeo_toc(&timer);
    stats.eval_grad_count++;
    return res;
}
double FatropAlg::eval_objective_curr()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    double res = 0.0;
    fatropnlp_->eval_obj(
        fatropdata_->obj_scale,
        fatropdata_->x_curr,
        res);
    stats.eval_obj_time += blasfeo_toc(&timer);
    stats.eval_obj_count++;
    return res;
}
fatrop_int FatropAlg::eval_dual_infeasiblity()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    fatropnlp_->eval_dual_inf(
        fatropdata_->obj_scale,
        fatropdata_->lam_curr,
        fatropdata_->grad_curr,
        fatropdata_->du_inf_curr);
    fatropdata_->eval_dual_inf_slack_eqs();
    stats.duinf_time += blasfeo_toc(&timer);
    return 0;
}
fatrop_int FatropAlg::perform_initializiation()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    fatrop_int res = fatropnlp_->initialize_slacks(
        fatropdata_->s_curr);
    res = fatropnlp_->initialize_dual(
        fatropdata_->grad_curr,
        fatropdata_->lam_calc,
        // fatropdata_->s_curr,
        fatropdata_->zL_curr,
        fatropdata_->zU_curr);
    stats.initialization_time += blasfeo_toc(&timer);
    return res;
}
fatrop_int FatropAlg::solve_pd_sys(double inertia_correction_w, double inertia_correction_c, double mu)
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    fatrop_int res = fatropnlp_->solve_pd_sys(
        inertia_correction_w,
        inertia_correction_c,
        fatropdata_->delta_x,
        fatropdata_->lam_calc,
        fatropdata_->delta_s,
        fatropdata_->sigma_total,
        fatropdata_->gradb_total);
    double el = blasfeo_toc(&timer);
    stats.compute_sd_time += el;
    return res;
}
