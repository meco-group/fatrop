/*
 * Fatrop - A fast trajectory optimization solver
 *  Copyright (C) 2022 - 2024 Lander Vanroye, KU Leuven. All rights reserved.
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
#include "fatrop/solver/FatropAlg.hpp"
using namespace fatrop;
using namespace std;

FatropAlg::FatropAlg(
    const shared_ptr<FatropNLP> &fatropnlp,
    const shared_ptr<FatropData> &fatropdata,
    const shared_ptr<FatropOptions> &fatropparams,
    const shared_ptr<Filter> &filter,
    const shared_ptr<LineSearch> &linesearch,
    const shared_ptr<Journaller> &journaller,
    const shared_ptr<FatropPrinter> &printer,
    const shared_ptr<FatropAlg> &orig_, const shared_ptr<FatropAlg> &resto_alg_, bool resto_problem)
    : fatropnlp_(fatropnlp),
      fatropdata_(fatropdata),
      fatropoptions_(fatropparams),
      filter_(filter),
      linesearch_(linesearch),
      journaller_(journaller),
      printer_(printer), orig_(orig_), resto_alg_(resto_alg_), resto_problem_(resto_problem)
{

    initialize();
    fatropnlp_->get_initial_sol_guess(fatropdata_->x_initial);
    fatropnlp->get_bounds(fatropdata->s_lower_orig, fatropdata->s_upper_orig);
}
void FatropAlg::initialize()
{
    // maxiter = fatropoptions_->max_iter;
    // kappa_d = fatropoptions_->kappa_d;
    fatropdata_->initialize();
    linesearch_->initialize();
}
void FatropAlg::reset()
{
    filter_->reset();
    if (!resto_problem_)
        fatropdata_->reset();
    fatropdata_->reset_caches();
    journaller_->reset();
    fatropnlp_->reset();
    linesearch_->reset();
    stats = FatropStats();
    if (!resto_problem_)
        start_iter_ = 0;
}
void FatropAlg::set_bounds(const vector<double> &lower, const vector<double> &upper)
{
    fatropdata_->set_bounds(lower, upper);
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
    fatropnlp_->update_mu(mu);
    double delta_w_last = 0.0;
    LineSearchInfo lsinfo;
    fatropnlp_->pre_solve(fatropdata_->x_curr, fatropdata_->s_curr);
    eval_constr_jac(); // todo twice evaluation
    eval_obj_grad_curr();
    if (!resto_problem_)
        fatropnlp_->initialize_slacks(mu0,
                                      fatropdata_->s_curr);
    if (!resto_problem_ && warm_start_init_point)
    {
        fatropdata_->warmstart_dual();
        fatropdata_->bound_z();
    }
    if (!resto_problem_ && !warm_start_init_point)
    {
        perform_initializiation_dual();
    }
    if (!resto_problem_)
        fatropdata_->bound_slacks();
    eval_constr_viol_curr();
    fatropdata_->theta_min = theta_min * MAX(1.0, fatropdata_->constr_viol_sum_curr());
    double theta_max = 1e4 * fatropdata_->constr_viol_sum_curr();
    filter_->augment(FilterData(0, std::numeric_limits<double>::infinity(), theta_max));
    fatrop_int ls = 0;
    double deltaw = 0;
    double deltac = 0.0;
    bool watch_dog_step = false;
    for (iter_count_ = start_iter_; iter_count_ < maxiter; iter_count_++)
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
        it_curr.iter = iter_count_;
        it_curr.mu = mu;
        it_curr.objective = eval_objective_curr();
        it_curr.constraint_violation = fatropdata_->constr_viol_max_curr();
        it_curr.du_inf = fatropdata_->dual_inf_max_curr();
        it_curr.ls = ls;
        it_curr.reg = deltaw;
        it_curr.resto = resto_problem_;
        if (no_no_full_steps_bc_filter >= 5)
        {
            bool reset_filter = (filter_reseted <= 5);
            if (reset_filter)
            {
                // printer_->level(1) << "resetted filter " << endl;
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
        journaller_->print_iterations(resto_problem_);
        double emu = fatropdata_->e_mu_curr(0.0);
        if (emu < acceptable_tol)
        {
            no_acceptable_steps++;
        }
        else
        {
            no_acceptable_steps = 0;
        }
        if (resto_problem_ && resto_stop_crit())
        {
            return 100;
        }

        if (emu < tol || (no_acceptable_steps >= acceptable_iter) || ((no_conse_small_sd == 2) && (mu <= mu_min)))
        // if (emu < tol)
        {
            double total_time = blasfeo_toc(&timer);
            journaller_->print_iterations(resto_problem_);
            if (no_conse_small_sd == 2)
            {
                printer_->level(1) << "WARNING fatrop returned bc of very small search direction" << endl;
            }
            if (emu > tol && no_acceptable_steps >= acceptable_iter)
            {
                printer_->level(1) << "WARNING fatrop returned acceptable tolerance" << endl;
            }
            stats.eval_cv_count += linesearch_->eval_cv_count;
            stats.eval_obj_count += linesearch_->eval_obj_count;
            stats.eval_cv_time += linesearch_->eval_cv_time;
            stats.eval_obj_time += linesearch_->eval_obj_time;
            stats.time_total = total_time;
            stats.iterations_count = iter_count_;
            stats.return_flag = 0;
            if (printer_->print_level() > 0)
            {
                stats.print(printer_->level(1));
            }
            if (resto_problem_)
                return 2;
            printer_->level(1) << "found solution" << endl;
            return 0;
        }
        // update mu
        // todo make a seperate class
        while (!watch_dog_step && mu > mu_min && (fatropdata_->e_mu_curr(mu) <= kappa_eta * mu || (no_conse_small_sd == 2)))
        {
            mu = MAX(mu_min, MIN(kappa_mu * mu, pow(mu, theta_mu)));
            fatropnlp_->update_mu(mu);
            if (resto_problem_)
                fatropdata_->obj_curr = eval_objective_curr();
            if (resto_problem_)
                eval_obj_grad_curr();
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
        // double stepsize = max(LinfScaled(fatropdata_->delta_x, fatropdata_->x_curr), LinfScaled(fatropdata_->delta_s, fatropdata_->s_curr));
        bool small_search_direction_curr = fatropdata_->small_step_size();
        lsinfo = linesearch_->find_acceptable_trial_point(mu, small_search_direction_curr || watch_dog_step, watch_dog_step);
        if (recalc_y && (deltac == 0.0) && (fatropdata_->constr_viol_max_curr() < recalc_y_feas_tol))
        {
            fatropnlp_->initialize_dual(
                fatropdata_->grad_curr_x,
                fatropdata_->grad_curr_s,
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
                // printer_->level(1) << "accepted watchdog step" << endl;
                watch_dog_step = false;
            }
            else
            {
                no_watch_dog_steps_taken++;
                if (no_watch_dog_steps_taken >= max_watchdog_steps)
                {
                    // reject watchdog step -- go back to x_k
                    // printer_->level(1) << "rejected watchdog step" << endl;
                    it_curr.type = 'x';
                    fatropdata_->restore_backup();
                    // delta_w_last = delta_w_last_backup;
                    watch_dog_step = false;
                    // todo make use of delta_x_backup and delta_s_backup
                    restore_watchdog_step = true;
                    iter_count_--;
                    continue;
                };
                it_curr.type = 'w';
            }
        }

        if (ls == 0)
        {
            // if already in resto phase: solver failed
            if (resto_problem_)
                return 1;
            // prepare for restoration phase
            int resto_res = start_resto_alg(mu, iter_count_);
            if (resto_res != 100)
            {
                stats.return_flag = 1;
                return 1;
            }
            else
            {
                return_from_resto_alg(mu);
            }
            iter_count_ = start_iter_ + 1;
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
    journaller_->print_iterations(resto_problem_);
    stats.return_flag = 1;
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
            fatropdata_->s_curr,
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
        fatropdata_->s_curr,
        fatropdata_->grad_curr_x,
        fatropdata_->grad_curr_s);
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
        fatropdata_->s_curr,
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
        fatropdata_->grad_curr_x,
        fatropdata_->grad_curr_s,
        fatropdata_->du_inf_curr, fatropdata_->du_inf_curr_s_wo_z);
    fatropdata_->eval_dual_inf_slack_eqs();
    stats.duinf_time += blasfeo_toc(&timer);
    return 0;
}
fatrop_int FatropAlg::perform_initializiation_dual()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    fatrop_int res = fatropnlp_->initialize_dual(
        fatropdata_->grad_curr_x,
        fatropdata_->grad_curr_s,
        fatropdata_->lam_calc,
        // fatropdata_->s_curr,
        fatropdata_->zL_curr,
        fatropdata_->zU_curr);
    if (res == 0 && fatropdata_->delta_dual_max() < lammax)
    {
        // printer_->level(1) << "accepted lam " << endl;
        fatropdata_->accept_dual_initializiaton();
    }
    else
    {
        // printer_->level(1) << "rejected lam " << endl;
        fatropdata_->lam_curr.SetConstant(0.0);
    }
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

fatrop_int FatropAlg::start_resto_alg(double mu, int iter)
{
    // augment filter with current iterate
    filter_->augment(FilterData(0, fatropdata_->obj_curr + fatropdata_->eval_barrier_func_curr(mu), fatropdata_->constr_viol_sum_curr()));
    fatrop_int n_ineqs = fatropdata_->n_ineqs;
    // set mu_init of resto alg
    mu = std::max(mu, fatropdata_->constr_viol_max_curr());
    resto_alg_->mu0 = mu;
    // set the starting iteration number
    resto_alg_->start_iter_ = iter + 1;
    // initialize primal variables
    resto_alg_->fatropdata_->x_curr.copy(fatropdata_->x_curr);
    // initialize the first part of the slack variables
    resto_alg_->fatropdata_->s_curr.block(0, fatropdata_->n_ineqs).copy(fatropdata_->s_curr);
    // set n and p variables to zero
    resto_alg_->fatropdata_->s_curr.block(fatropdata_->n_ineqs, 2 * fatropdata_->n_ineqs) = 0.;
    // evaluate the constraint jacobian of resto_alg
    // the constraint jacobian is required for the slack initialization (todo: make this more efficient)
    resto_alg_->eval_constr_jac();
    // call initialize slecks from the resto nlp
    resto_alg_->fatropnlp_->initialize_slacks(mu, resto_alg_->fatropdata_->s_curr);
    // initialize equality multipliers
    resto_alg_->fatropdata_->lam_curr = 0.;
    // initialize z0, zn and zp
    for (fatrop_int i = 0; i < n_ineqs; i++)
    {
        resto_alg_->fatropdata_->zL_curr.at(i) = std::min(1000., fatropdata_->zL_curr.at(i));
        resto_alg_->fatropdata_->zU_curr.at(i) = std::min(1000., fatropdata_->zU_curr.at(i));
        resto_alg_->fatropdata_->zL_curr.at(n_ineqs + i) = mu / resto_alg_->fatropdata_->s_curr.at(n_ineqs + i);
        resto_alg_->fatropdata_->zL_curr.at(2 * n_ineqs + i) = mu / resto_alg_->fatropdata_->s_curr.at(2 * n_ineqs + i);
    }
    // call resto alg
    return resto_alg_->optimize();
}

void FatropAlg::update_options(const FatropOptions &options)
{
    maxiter = options.max_iter.get();
    kappa_d = options.kappa_d.get();
    tol = options.tol.get();
    acceptable_tol = options.acceptable_tol.get();
    max_watchdog_steps = options.max_watchdog_steps.get();
    acceptable_iter = options.acceptable_iter.get();
    lammax = options.lammax.get();
    mu0 = options.mu_init.get();
    kappa_eta = options.kappa_eta.get();
    kappa_mu = options.kappa_mu.get();
    theta_mu = options.theta_mu.get();
    delta_w0 = options.delta_w0.get();
    delta_wmin = options.delta_wmin.get();
    kappa_wmin = options.kappa_wmin.get();
    kappa_wplus = options.kappa_wplus.get();
    kappa_wplusem = options.kappa_wplusem.get();
    delta_c_stripe = options.delta_c_stripe.get();
    kappa_c = options.kappa_c.get();
    warm_start_init_point = options.warm_start_init_point.get();
    theta_min = options.theta_min.get();
    recalc_y = options.recalc_y.get();
    recalc_y_feas_tol = options.recalc_y_feas_tol.get();
};

fatrop_int FatropAlg::return_from_resto_alg(double mu)
{
    // compute delta x
    axpby(1.0, resto_alg_->fatropdata_->x_curr, -1.0, fatropdata_->x_curr, fatropdata_->delta_x);
    // compute delta s
    axpby(1.0, resto_alg_->fatropdata_->s_curr.block(0, fatropdata_->n_ineqs), -1.0, fatropdata_->s_curr, fatropdata_->delta_s);
    // compute delta z
    fatropdata_->compute_delta_z();
    // compute maximum step size
    double alpha_primal = 1.0;
    double alpha_dual = 1.0;
    fatropdata_->maximum_step_size(alpha_primal, alpha_dual, std::max(1 - mu, 0.99));
    // update trial step
    fatropdata_->update_trial_step(alpha_primal, alpha_dual);
    // augment filter
    // filter_->augment(FilterData(0, linesearch_->eval_obj_trial() + fatropdata_->eval_barrier_func_trial(mu), fatropdata_->constr_viol_sum_curr()));
    // evaluate some quantities before accepting the trial step
    linesearch_->eval_constr_viol_trial();
    // accept trial step
    fatropdata_->accept_trial_step();
    // initialize equality multipliers
    eval_constr_jac(); // todo twice evaluation
    eval_obj_grad_curr();
    perform_initializiation_dual();
    fatropdata_->modify_dual_bounds(mu);
    // update iteration number
    start_iter_ = resto_alg_->iter_count_ - 1;
    return 0;
}
bool FatropAlg::resto_stop_crit()
{
    auto orig_p = orig_.lock();
    // evaluate original's problem constraint violation of the current solution
    orig_p->fatropnlp_->eval_constraint_viol(fatropdata_->x_curr, fatropdata_->s_curr, orig_p->fatropdata_->g_next);
    double cv_i = orig_p->fatropdata_->constr_viol_sum_next();
    double obj_i = 0.0;
    orig_p->fatropnlp_->eval_obj(fatropdata_->obj_scale, fatropdata_->x_curr, fatropdata_->s_curr, obj_i);
    // check if acceptable to original filter
    return orig_p->filter_->is_acceptable(FilterData(0, obj_i, cv_i)) && cv_i < 0.9 * orig_p->fatropdata_->constr_viol_sum_curr();
}