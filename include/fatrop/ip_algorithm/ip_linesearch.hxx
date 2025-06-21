//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_linesearch_hxx__
#define __fatrop_ip_algorithm_ip_linesearch_hxx__
#include "fatrop/common/options.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_linesearch.hpp"
#include "fatrop/ip_algorithm/ip_resto_phase.hpp"
#include "fatrop/ip_algorithm/ip_utils.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include <cmath>
#include <limits>

namespace fatrop
{
    template <typename LinearSolverType, typename ProblemType>
    IpLinesearch<LinearSolverType, ProblemType>::IpLinesearch(
        const IpDataSp &ipdata, const LinearSolverSp &linear_solver,
        const RestorationPhaseSp &restoration_phase)
        : ipdata_(ipdata), linear_solver_(linear_solver), restoration_phase_(restoration_phase),
          soc_rhs_x_(ipdata->current_iterate().nlp()->nlp_dims().number_of_variables),
          soc_rhs_s_(ipdata->current_iterate().nlp()->nlp_dims().number_of_ineq_constraints),
          soc_rhs_g_(ipdata->current_iterate().nlp()->nlp_dims().number_of_eq_constraints),
          soc_g_accumulate_(ipdata->current_iterate().nlp()->nlp_dims().number_of_eq_constraints),
          soc_rhs_cl_(ipdata->current_iterate().nlp()->nlp_dims().number_of_ineq_constraints),
          soc_rhs_cu_(ipdata->current_iterate().nlp()->nlp_dims().number_of_ineq_constraints)
    {
        // avoid dynamic memory allocation during solving
        filter().reserve(max_iter_ + 1);
    }

    template <typename LinearSolverType, typename ProblemType>
    void IpLinesearch<LinearSolverType, ProblemType>::init_this_line_search(bool in_watchdog)
    {
        IpIterateType &curr_it = ipdata_->current_iterate();
        if (!in_watchdog)
        {
            reference_theta_ = norm_l1(curr_it.constr_viol());
            reference_barr_ = curr_it.barrier_value() + curr_it.obj_value();
            reference_grad_bar_delta_ =
                curr_it.linear_decrease_barrier() + curr_it.linear_decrease_objective();
        }
        else
        {
            reference_theta_ = watchdog_theta_;
            reference_barr_ = watchdog_barr_;
            reference_grad_bar_delta_ = watchdog_grad_bar_delta_;
        }
    }
    template <typename LinearSolverType, typename ProblemType>
    bool IpLinesearch<LinearSolverType, ProblemType>::do_backtracking_line_search(
        bool skip_first_trial_point, Scalar &alpha_primal, bool &corr_taken, bool &soc_taken,
        Index &n_steps, bool &evaluation_error)
    {
        IpIterateType &curr_it = ipdata_->current_iterate();
        evaluation_error = false;
        bool accept = false;
        Scalar alpha_primal_max = curr_it.maximum_step_size_primal(curr_it.tau());
        Scalar alpha_min = in_watchdog_ ? alpha_primal_max : compute_alpha_min();
        // start the line search with the maximum step size
        alpha_primal = alpha_primal_max;
        // alpha used for ftype and amijo test
        const Scalar alpha_primal_test = in_watchdog_ ? watchdog_alpha_primal_test_ : alpha_primal;
        if (skip_first_trial_point)
        {
            alpha_primal *= alpha_red_factor_;
        }
        // todo corrector step not implemented
        while (alpha_primal > alpha_min || n_steps == 0)
        {
            // set the trial point
            IpIterateType &trial_it = ipdata_->trial_iterate();
            trial_it.set_primal_x(curr_it.primal_x() + alpha_primal * curr_it.delta_primal_x());
            trial_it.set_primal_s(curr_it.primal_s() + alpha_primal * curr_it.delta_primal_s());
            bool eval_success =
                (trial_it.constr_viol().is_finite() && std::isfinite(trial_it.obj_value()) &&
                 std::isfinite(trial_it.barrier_value()));
            // todo not implemented: accept every trial step, accept after max steps
            // check if we can evaluate the objective, barrier and constraints
            if (eval_success)
                accept = check_acceptability_of_trial_point(alpha_primal_test);
            if (!eval_success)
            {
                evaluation_error = true;
                accept = false;
            }
            if (accept || in_watchdog_)
                break;
            // second order correction
            if (!evaluation_error && alpha_primal == alpha_primal_max)
            {
                Scalar theta_curr = norm_l1(curr_it.constr_viol());
                Scalar theta_trial = norm_l1(trial_it.constr_viol());
                if (theta_curr <= theta_trial)
                {
                    accept = try_second_order_correction(alpha_primal_test, alpha_primal);
                }
                if (accept)
                {
                    soc_taken = true;
                    break;
                }
            }
            alpha_primal *= alpha_red_factor_;
            n_steps++;
        }

        Scalar alpha_dual_max = curr_it.maximum_step_size_dual(curr_it.tau());

        char info_alpha_primal_char = '?';
        if (!accept && in_watchdog_)
        {
            info_alpha_primal_char = 'w';
        }
        else if (accept)
        {
            info_alpha_primal_char = update_for_next_iteration(alpha_primal_max);
        }
        if (soc_taken)
        {
            info_alpha_primal_char = (char)toupper(info_alpha_primal_char);
        }
        update_step_info(alpha_primal, alpha_dual_max, n_steps, info_alpha_primal_char);
        return accept;
    }

    template <typename LinearSolverType, typename ProblemType>
    void IpLinesearch<LinearSolverType, ProblemType>::reset()
    {
        theta_min_ = -1.;
        theta_max_ = -1.;
        filter_reset_count_ = 0;
        in_watchdog_ = false;
        tiny_step_last_iteration_ = false;
        sucessive_tiny_step_count_ = 0;
        acceptable_iteration_count_ = -1;
        last_mu_ = -1.;
        reset_linesearch();
        watchdog_shortened_iter_count_ = 0;
        watchdog_trial_iter_count_ = 0;
        in_soft_resto_phase_ = false;
        soft_resto_counter_ = 0;
        if (restoration_phase_)
            restoration_phase_->reset();
    }

    /**
     * this is a seperate function because not all quanities have to be resetted when called from mu
     * updater
     */

    template <typename LinearSolverType, typename ProblemType>
    void IpLinesearch<LinearSolverType, ProblemType>::reset_linesearch()
    {
        filter_reject_count_ = 0;
        last_rejection_due_to_filter_ = false;
        filter().reset();
    }

    template <typename LinearSolverType, typename ProblemType>
    IpFilter &IpLinesearch<LinearSolverType, ProblemType>::filter()
    {
        return filter_;
    }

    template <typename LinearSolverType, typename ProblemType>
    bool IpLinesearch<LinearSolverType, ProblemType>::find_acceptable_trial_point()
    {
        IpIterateType *curr_it = &ipdata_->current_iterate();
        IpIterateType *trial_it = &ipdata_->trial_iterate();
        Scalar curr_mu = (*curr_it).mu();
        bool mu_changed = (curr_mu != last_mu_);
        if (mu_changed) // or first call
        {
            // inactivate the watchdog
            in_watchdog_ = false;
            watchdog_shortened_iter_count_ = 0;
        }
        last_mu_ = curr_mu;
        // todo: ipopt has a kind of fallback mechanism which falls back to an acceptable iterate,
        // we do not have implemented this yet
        // todo: we don´t have emergency mode - fall back mechanism implmented
        bool goto_resto = false;
        init_this_line_search(in_watchdog_);
        // todo "start with resto" not implemented
        // todo "expect infeasible problem" not implemented
        bool accept = false;
        Index n_steps = 0;
        Scalar alpha_primal = 0.;
        bool tiny_step = (!goto_resto && detect_tiny_step());
        if (in_watchdog_ && (goto_resto || tiny_step))
        {
            stop_watchdog();
            curr_it = &ipdata_->current_iterate();
            trial_it = &ipdata_->trial_iterate();
            goto_resto = false;
            tiny_step = false;
        }
        // Check if we want to wake up the watchdog
        if (watchdog_shortened_iter_trigger_ > 0 && !in_watchdog_ && !goto_resto && !tiny_step &&
            !in_soft_resto_phase_ /*&& !expect_infeasible_problem_ */ &&
            watchdog_shortened_iter_count_ >= watchdog_shortened_iter_trigger_)
        {
            start_watchdog();
            curr_it = &ipdata_->current_iterate();
            trial_it = &ipdata_->trial_iterate();
        }
        //  handle the case of tiny steps
        if (tiny_step)
        {
            Scalar alpha_primal = (*curr_it).maximum_step_size_primal((*curr_it).tau());
            (*trial_it).set_primal_x((*curr_it).primal_x() +
                                     alpha_primal * (*curr_it).delta_primal_x()); // x
            (*trial_it).set_primal_s((*curr_it).primal_s() +
                                     alpha_primal * (*curr_it).delta_primal_s()); // s
            // let's see if we can evaluate objective barrier and contraints and if they are finite
            if (!(*trial_it).constr_viol().is_finite() || !std::isfinite((*trial_it).obj_value()) ||
                !std::isfinite((*trial_it).barrier_value()))
            {
                // dont accept the tiny step
                tiny_step = false;
            }
            if (tiny_step && tiny_step_last_iteration_)
                ipdata_->set_tiny_step_flag(true);

            Scalar delta_y_norm = norm_inf((*trial_it).delta_dual_eq());
            if (delta_y_norm < tiny_step_y_tol_)
            {
                tiny_step_last_iteration_ = false;
            }
            else
            {
                tiny_step_last_iteration_ = true;
            }

            accept = true;
        }
        else
        {
            tiny_step_last_iteration_ = false;
        }
        if (!goto_resto && !tiny_step)
        {
            if (in_soft_resto_phase_)
            {
                soft_resto_counter_++;
                if (soft_resto_counter_ > max_soft_resto_iters_)
                {
                    accept = false;
                }
                else
                {
                    bool satisfies_original_criterion = false;
                    accept = try_soft_resto_step(satisfies_original_criterion);
                    if (accept)
                    {
                        (*trial_it).step_info().alpha_primal_type = 's';
                        if (satisfies_original_criterion)
                        {
                            in_soft_resto_phase_ = false;
                            soft_resto_counter_ = 0;
                            (*trial_it).step_info().alpha_primal_type = 'S';
                        }
                    }
                }
            }
            else
            {
                bool done = false;
                bool skip_first_trial_point = false;
                bool evaluation_error = false;
                while (!done)
                {
                    bool corr_taken = false;
                    bool soc_taken = false;
                    accept = do_backtracking_line_search(skip_first_trial_point, alpha_primal,
                                                         corr_taken, soc_taken, n_steps,
                                                         evaluation_error);
                    if (in_watchdog_ && accept)
                    {
                        in_watchdog_ = false;
                        done = true;
                    }
                    if (in_watchdog_ && !accept)
                    {
                        watchdog_trial_iter_count_++;
                        if (evaluation_error ||
                            watchdog_trial_iter_count_ >= watchdog_trial_iter_max_)
                        {
                            if (in_watchdog_)
                            {
                                stop_watchdog();
                                curr_it = &ipdata_->current_iterate();
                                trial_it = &ipdata_->trial_iterate();
                            }
                            skip_first_trial_point = true;
                        }
                        else
                        {
                            done = true;
                            accept = true;
                        }
                    }
                    if (!in_watchdog_)
                    {
                        done = true;
                    }
                }
            } // else if(in_soft_restoration_phase)
        } // if(!goto_resto && !tiny_step)
        // if the line search is unsuccesfull go to the restoration phase
        if (!accept)
        {
            // todo rigorous mode not implemented
            if (!in_soft_resto_phase_ && !goto_resto /* && ! expect_infeasible_problem*/)
            {
                prepare_resto_phase_start();
                bool satisfies_original_criterion = false;
                accept = try_soft_resto_step(satisfies_original_criterion);
                if (accept)
                {
                    if (satisfies_original_criterion)
                    {
                        (*trial_it).step_info().alpha_primal_type = 'S';
                    }
                    else
                    {
                        (*trial_it).step_info().alpha_primal_type = 's';
                        in_soft_resto_phase_ = true;
                    }
                }
            }
            if (!accept) // go to restoration phase
            {
                // f_out << "calling restoration phase" << std::endl;
                fatrop_assert_msg(restoration_phase_,
                                  "Restoration phase not set in line search object. Maybe called "
                                  "from restoration phase algorithm.");
                /**
                 * todo set the step info
                 */
                accept = restoration_phase_->perform_restoration();
                PRINT_DEBUG << "Restoration phase failed" << std::endl;

                in_soft_resto_phase_ = false;
                soft_resto_counter_ = 0;
                watchdog_shortened_iter_count_ = 0;
            }
        }
        else if (!in_soft_resto_phase_ || tiny_step)
        {
            // no restoration phase and now updating dual variables
            if (accept)
            {
                auto alpha_max = (*curr_it).maximum_step_size((*curr_it).tau());
                Scalar alpha_primal_max = alpha_max.first;
                Scalar alpha_dual_max = alpha_max.second;
                perform_dual_step(alpha_primal, alpha_dual_max);

                if (n_steps == 0)
                {
                    watchdog_shortened_iter_count_ = 0;
                }
                if (n_steps > 0)
                {
                    watchdog_shortened_iter_count_++;
                }
            }
        }
        return accept;
    }
    template <typename LinearSolverType, typename ProblemType>
    void IpLinesearch<LinearSolverType, ProblemType>::update_step_info(
        const Scalar alpha_primal, const Scalar alpha_dual, const Index n_steps,
        const char info_alpha_primal_char)
    {
        IpIterateType &curr_it = ipdata_->current_iterate();
        IpIterateType &trial_it = ipdata_->trial_iterate();

        trial_it.step_info().alpha_primal = alpha_primal;
        trial_it.step_info().alpha_dual = alpha_dual;
        trial_it.step_info().step_length =
            std::max(norm_inf(trial_it.delta_primal_x()), norm_inf(trial_it.delta_primal_s()));
        trial_it.step_info().inertia_correction_primal =
            curr_it.search_dir_info().inertia_correction_primal;
        trial_it.step_info().inertia_correction_dual =
            curr_it.search_dir_info().inertia_correction_dual;
        trial_it.step_info().ls_iter = n_steps;
        trial_it.step_info().alpha_primal_type = info_alpha_primal_char;
    }

    template <typename LinearSolverType, typename ProblemType>
    void IpLinesearch<LinearSolverType, ProblemType>::prepare_resto_phase_start()
    {
        // augment the filter with the current iterate
        augment_filter();
    }
    template <typename LinearSolverType, typename ProblemType>
    char IpLinesearch<LinearSolverType, ProblemType>::update_for_next_iteration(
        const Scalar alpha_primal_test)
    {
        char info_alpha_primal_char;
        if (!is_f_type(alpha_primal_test) || !armijo_holds(alpha_primal_test))
        {
            info_alpha_primal_char = 'h';
        }
        else
        {
            info_alpha_primal_char = 'f';
        }
        return info_alpha_primal_char;
    }

    template <typename LinearSolverType, typename ProblemType>
    void IpLinesearch<LinearSolverType, ProblemType>::start_watchdog()
    {
        fatrop_assert_msg(!in_watchdog_, "start watchdog called but watchdog is already active");
        in_watchdog_ = true;
        watchdog_trial_iter_count_ = 0;
        IpIterateType &curr_it = ipdata_->current_iterate();
        watchdog_theta_ = norm_l1(curr_it.constr_viol());
        watchdog_barr_ = curr_it.barrier_value() + curr_it.obj_value();
        watchdog_grad_bar_delta_ =
            curr_it.linear_decrease_barrier() + curr_it.linear_decrease_objective();
        watchdog_alpha_primal_test_ = curr_it.maximum_step_size(curr_it.tau()).first;
        ipdata_->store_current_iterate();
    }

    template <typename LinearSolverType, typename ProblemType>
    void IpLinesearch<LinearSolverType, ProblemType>::stop_watchdog()
    {
        fatrop_assert_msg(in_watchdog_, "stop watchdog called but watchdog is not active");
        in_watchdog_ = false;
        reference_theta_ = watchdog_theta_;
        reference_barr_ = watchdog_barr_;
        reference_grad_bar_delta_ = watchdog_grad_bar_delta_;
        ipdata_->restore_current_iterate();
    }

    template <typename LinearSolverType, typename ProblemType>
    void IpLinesearch<LinearSolverType, ProblemType>::perform_dual_step(const Scalar alpha_primal,
                                                                        const Scalar alpha_dual)
    {
        IpIterateType& curr_it = ipdata_->current_iterate();
        IpIterateType& trial_it = ipdata_->trial_iterate();
        trial_it.set_dual_eq(curr_it.dual_eq() + alpha_primal * curr_it.delta_dual_eq());
        trial_it.set_dual_bounds_l(curr_it.dual_bounds_l() +
                                   alpha_dual * curr_it.delta_dual_bounds_l());
        trial_it.set_dual_bounds_u(curr_it.dual_bounds_u() +
                                   alpha_dual * curr_it.delta_dual_bounds_u());
        trial_it.modify_dual_bounds(curr_it.mu());
    }

    template <typename LinearSolverType, typename ProblemType>
    bool IpLinesearch<LinearSolverType, ProblemType>::try_second_order_correction(
        const Scalar alpha_primal_test, Scalar &alpha_primal)
    {
        if (max_soc_ == 0)
            return false;
        IpIterateType &curr_it = ipdata_->current_iterate();
        IpIterateType &trial_it = ipdata_->trial_iterate();

        bool accept = false;
        Index count_soc = 0;
        Scalar alpha_primal_soc = alpha_primal;
        Scalar theta_soc_old = 0.;
        Scalar theta_trial = norm_l1(trial_it.constr_viol());
        soc_g_accumulate_ = curr_it.constr_viol();

        while (count_soc < max_soc_ && !accept &&
               (count_soc == 0 || theta_trial < kappa_soc_ * theta_soc_old))
        {
            theta_soc_old = theta_trial;

            // set up the rhs for the second order correction
            soc_g_accumulate_ = soc_g_accumulate_ + alpha_primal_soc * trial_it.constr_viol();

            soc_rhs_x_ = curr_it.dual_infeasibility_x();
            soc_rhs_s_ = curr_it.dual_infeasibility_s();
            soc_rhs_g_ = soc_g_accumulate_;
            Scalar mu = curr_it.mu();
            soc_rhs_cl_.block(soc_rhs_cl_.m(), 0) = curr_it.relaxed_complementarity_l();
            soc_rhs_cu_.block(soc_rhs_cu_.m(), 0) = curr_it.relaxed_complementarity_u();

            // solve the linear system
            typename LinearSolverType::LinearSystemType ls(
                curr_it.info(), curr_it.jacobian(), curr_it.hessian(), curr_it.Dx(),
                curr_it.De_is_zero(), curr_it.De(), curr_it.delta_lower(), curr_it.delta_upper(),
                curr_it.dual_bounds_l(), curr_it.dual_bounds_u(), soc_rhs_x_, soc_rhs_s_,
                soc_rhs_g_, soc_rhs_cl_, soc_rhs_cu_);

            LinsolReturnFlag ret = linear_solver_->solve_in_place_rhs(ls);
            // check if the linear system was solved successfully
            if (!(ret == LinsolReturnFlag::SUCCESS || ret == LinsolReturnFlag::ITREF_MAX_ITER ||
                  ret == LinsolReturnFlag::ITREF_INCREASE))
                break;
            // compute the maximum step size for the second order correction
            alpha_primal_soc = curr_it.maximum_step_size_primal(curr_it.tau(), soc_rhs_s_);
            // set the new trial iterate
            trial_it.set_primal_x(curr_it.primal_x() + alpha_primal_soc * soc_rhs_x_);
            trial_it.set_primal_s(curr_it.primal_s() + alpha_primal_soc * soc_rhs_s_);
            // check if the trial point is acceptable
            accept = check_acceptability_of_trial_point(alpha_primal_test);
            if (accept)
            {
                // set the new search directions
                curr_it.set_delta_primal_x(soc_rhs_x_);
                curr_it.set_delta_primal_s(soc_rhs_s_);
                curr_it.set_delta_dual_eq(soc_rhs_g_);
                curr_it.set_delta_dual_bounds_l(soc_rhs_cl_);
                curr_it.set_delta_dual_bounds_u(soc_rhs_cu_);
                alpha_primal = alpha_primal_soc;
            }
            else
            {
                theta_trial = norm_l1(trial_it.constr_viol());
                count_soc++;
            }
        }
        return accept;
    }
    // unnamed namespace for limiting scope to compilation unit of this file
    namespace
    {
        template <typename ProblemType> Scalar primal_dual_error(IpIterate<ProblemType> &iterate)
        {
            Scalar res = 0.;
            res += norm_l1_divided_by_m(iterate.dual_infeasibility_x());
            res += norm_l1_divided_by_m(iterate.dual_infeasibility_s());
            res += norm_l1_divided_by_m(iterate.constr_viol());
            res += norm_l1_divided_by_m(iterate.relaxed_complementarity_l());
            res += norm_l1_divided_by_m(iterate.relaxed_complementarity_u());
            return res;
        }
    }
    template <typename LinearSolverType, typename ProblemType>
    bool IpLinesearch<LinearSolverType, ProblemType>::try_soft_resto_step(
        bool &satisfies_original_criterion)
    {
        if (soft_rest_pd_error_reduction_factor_ == 0.)
            return false;

        satisfies_original_criterion = false;
        IpIterateType &curr_it = ipdata_->current_iterate();
        IpIterateType &trial = ipdata_->trial_iterate();
        Scalar alpha_primal_max = curr_it.maximum_step_size_primal(curr_it.tau());
        Scalar alpha_dual_max = curr_it.maximum_step_size_dual(curr_it.tau());
        Scalar alpha = std::min(alpha_primal_max, alpha_dual_max);
        // todo: in contrast to Ipopt, we dont check for evaluation errors, we will however throw an
        // exception at check_acceptability_of_trial_point but we dont have a mechanism to catch
        // this exception
        trial.set_primal_x(curr_it.primal_x() + alpha * curr_it.delta_primal_x());
        trial.set_primal_s(curr_it.primal_s() + alpha * curr_it.delta_primal_s());
        perform_dual_step(alpha, alpha);
        if (check_acceptability_of_trial_point(0.))
        {
            satisfies_original_criterion = true;
            return true;
        }

        Scalar mu = curr_it.mu();
        Scalar trial_pd_error = primal_dual_error(trial);
        Scalar curr_pd_error = primal_dual_error(curr_it);
        // check if there is sufficient reduction in the primal dual error
        if (trial_pd_error <= soft_rest_pd_error_reduction_factor_ * curr_pd_error)
        {
            update_step_info(alpha, alpha, 0, 's');
            return true;
        }
        return false;
    }

    template <typename LinearSolverType, typename ProblemType>
    bool IpLinesearch<LinearSolverType, ProblemType>::is_acceptable_to_current_iterate(
        const Scalar trial_barr, const Scalar trial_theta, const bool called_from_resto /*=false*/)
    {
        using namespace internal;
        // check if the barrier objective function is increasing too rapidly
        if (!called_from_resto && trial_barr > reference_barr_)
        {
            Scalar bas_val = 1.;
            if (std::abs(reference_barr_) > 10.)
                bas_val = std::log10(std::abs(reference_barr_));
            if (std::log10(trial_barr - reference_barr_) > obj_max_incr_ + bas_val)
                return false;
        }
        return (compare_le(trial_theta, (1. - gamma_theta_) * reference_theta_, reference_theta_) ||
                compare_le(trial_barr - reference_barr_, -gamma_phi_ * reference_theta_,
                           reference_barr_));
    }

    template <typename LinearSolverType, typename ProblemType>
    bool
    IpLinesearch<LinearSolverType, ProblemType>::is_acceptable_to_filter(const Scalar trial_barr,
                                                                         const Scalar trial_theta)
    {
        return filter().is_acceptable({trial_barr, trial_theta}); // Placeholder return
    }

    template <typename LinearSolverType, typename ProblemType>
    bool IpLinesearch<LinearSolverType, ProblemType>::armijo_holds(const Scalar alpha_primal)
    {
        using namespace internal;
        Scalar trial_barr =
            ipdata_->trial_iterate().barrier_value() + ipdata_->trial_iterate().obj_value();
        return compare_le(trial_barr - reference_barr_,
                          eta_phi_ * alpha_primal * reference_grad_bar_delta_, reference_barr_);
    }

    template <typename LinearSolverType, typename ProblemType>
    bool IpLinesearch<LinearSolverType, ProblemType>::detect_tiny_step()
    {
        if (tiny_step_tol_ == 0.)
            return false;
        const VecRealView &curr_x = ipdata_->current_iterate().primal_x();
        const VecRealView &delta_x = ipdata_->current_iterate().delta_primal_x();
        const Scalar max_step_x = max_el(abs(delta_x) / (abs(curr_x) + 1.));
        if (max_step_x > tiny_step_tol_)
            return false;
        const VecRealView &curr_s = ipdata_->current_iterate().primal_s();
        const VecRealView &delta_s = ipdata_->current_iterate().delta_primal_s();
        const Scalar max_step_s = max_el(abs(delta_s) / (abs(curr_s) + 1.));
        if (max_step_s > tiny_step_tol_)
            return false;
        const Scalar cviol = norm_l1(ipdata_->current_iterate().constr_viol());
        if (cviol > 1e-4)
            return false;
        return true;
    }

    template <typename LinearSolverType, typename ProblemType>
    bool IpLinesearch<LinearSolverType, ProblemType>::is_f_type(const Scalar alpha_primal)
    {
        Scalar mach_eps = std::numeric_limits<Scalar>::epsilon();
        if (reference_theta_ == 0. && reference_grad_bar_delta_ > 0. &&
            reference_grad_bar_delta_ < 100. * mach_eps)
        {
            reference_grad_bar_delta_ = -mach_eps;
        }
        return (reference_grad_bar_delta_ < 0. &&
                alpha_primal * std::pow(-reference_grad_bar_delta_, s_phi_) >
                    delta_ * std::pow(reference_theta_, s_theta_));
    }

    template <typename LinearSolverType, typename ProblemType>
    void IpLinesearch<LinearSolverType, ProblemType>::augment_filter()
    {
        Scalar phi_add = reference_barr_ - gamma_phi_ * reference_theta_;
        Scalar theta_add = (1. - gamma_theta_) * reference_theta_;
        filter().augment({phi_add, theta_add});
    }
    template <typename LinearSolverType, typename ProblemType>
    Scalar IpLinesearch<LinearSolverType, ProblemType>::compute_alpha_min()
    {
        IpIterateType &curr_it = ipdata_->current_iterate();
        Scalar gBD = curr_it.linear_decrease_barrier() + curr_it.linear_decrease_objective();
        Scalar curr_theta = norm_l1(curr_it.constr_viol());
        Scalar alpha_min = gamma_theta_;
        if (gBD < 0.)
        {
            alpha_min = std::min(gamma_theta_, gamma_phi_ * curr_theta / (-gBD));
            if (curr_theta <= theta_min_)
            {
                alpha_min = std::min(alpha_min, delta_ * std::pow(curr_theta, s_theta_) /
                                                    std::pow(-gBD, s_phi_));
            }
        }
        return alpha_min_frac_ * alpha_min;
    }

    template <typename LinearSolverType, typename ProblemType>
    bool IpLinesearch<LinearSolverType, ProblemType>::check_acceptability_of_trial_point(
        const Scalar alpha_primal)
    {
        bool accept;
        IpIterateType &trial_it = ipdata_->trial_iterate();
        Scalar trial_theta = norm_l1(trial_it.constr_viol());
        if (theta_max_ < 0.)
        {
            theta_max_ = theta_max_fact_ * std::max(1.0, reference_theta_);
        }
        if (theta_min_ < 0.)
        {
            theta_min_ = theta_min_fact_ * std::max(1.0, reference_theta_);
        }
        // if constraint violation is too large -> reject
        if (theta_max_ > 0 && trial_theta > theta_max_)
            return false;
        Scalar trial_barr = trial_it.barrier_value() + trial_it.obj_value();
        fatrop_assert_msg(std::isfinite(trial_barr), "trial_barr is not finite");
        // check if a poin tis acceptable wrt current iterate
        if (alpha_primal > 0. && is_f_type(alpha_primal) && reference_theta_ < theta_min_)
            accept = armijo_holds(alpha_primal);
        else
            accept = is_acceptable_to_current_iterate(trial_barr, trial_theta);

        last_rejection_due_to_filter_ = false;
        if (!accept)
        {
            return false;
        }
        // now check if the point is acceptable wrt the filter
        accept = is_acceptable_to_filter(trial_barr, trial_theta);
        if (!accept)
        {
            last_rejection_due_to_filter_ = true;
            return false;
        }
        // accept == true
        // Filter reset heuristic
        if (max_filter_resets_ <= 0 || filter_reset_count_ <= max_filter_resets_ ||
            filter_reject_count_ <= filter_reset_trigger_)
        {
            return true;
        }

        if (!last_rejection_due_to_filter_)
        {
            filter_reject_count_ = 0;
            return true;
        }
        else
        {
            filter_reject_count_++;
        }

        if (filter_reject_count_ >= filter_reset_trigger_)
        {
            filter().reset();
            // filter_reset_count_++; // TODO: Ipopt implementation does not seem to
            // increment this while I would expect it to
        }
        return true;
    }

    template <typename LinearSolverType, typename ProblemType>
    void IpLinesearch<LinearSolverType, ProblemType>::register_options(OptionRegistry &registry)
    {
        registry.register_option("max_soc", &IpLinesearch::set_max_soc, this);
        registry.register_option("kappa_soc", &IpLinesearch::set_kappa_soc, this);
        registry.register_option("theta_min", &IpLinesearch::set_theta_min, this);
        registry.register_option("theta_max", &IpLinesearch::set_theta_max, this);
        registry.register_option("theta_min_fact", &IpLinesearch::set_theta_min_fact, this);
        registry.register_option("s_phi", &IpLinesearch::set_s_phi, this);
        registry.register_option("s_theta", &IpLinesearch::set_s_theta, this);
        registry.register_option("delta", &IpLinesearch::set_delta, this);
        registry.register_option("tiny_step_tol", &IpLinesearch::set_tiny_step_tol, this);
        registry.register_option("tiny_step_y_tol", &IpLinesearch::set_tiny_step_y_tol, this);
        registry.register_option("gamma_theta", &IpLinesearch::set_gamma_theta, this);
        registry.register_option("gamma_phi", &IpLinesearch::set_gamma_phi, this);
        registry.register_option("eta_phi", &IpLinesearch::set_eta_phi, this);
        registry.register_option("alpha_min_frac", &IpLinesearch::set_alpha_min_frac, this);
        registry.register_option("theta_max_fact", &IpLinesearch::set_theta_max_fact, this);
        registry.register_option("max_filter_resets", &IpLinesearch::set_max_filter_resets, this);
        registry.register_option("filter_reset_trigger", &IpLinesearch::set_filter_reset_trigger,
                                 this);
        registry.register_option("obj_max_incr", &IpLinesearch::set_obj_max_incr, this);
        registry.register_option("watchdog_shortened_iter_trigger",
                                 &IpLinesearch::set_watchdog_shortened_iter_trigger, this);
        registry.register_option("watchdog_trial_iter_max",
                                 &IpLinesearch::set_watchdog_trial_iter_max, this);
        registry.register_option("alpha_red_factor", &IpLinesearch::set_alpha_red_factor, this);
        registry.register_option("max_iter", &IpLinesearch::set_max_iter, this);
        registry.register_option("soft_rest_pd_error_reduction_factor",
                                 &IpLinesearch::set_soft_rest_pd_error_reduction_factor, this);
        registry.register_option("max_soft_resto_iters", &IpLinesearch::set_max_soft_resto_iters,
                                 this);
    }

} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_linesearch_hxx__
