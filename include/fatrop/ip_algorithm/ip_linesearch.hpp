//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_linesearch_hpp__
#define __fatrop_ip_algorithm_ip_linesearch_hpp__
#include "fatrop/common/fwd.hpp"
#include "fatrop/context/context.hpp"
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/ip_algorithm/ip_filter.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/nlp/fwd.hpp"
#include <limits>
#include <memory>

namespace fatrop
{
    /**
     * @brief Base class for line search algorithms in interior point methods.
     */
    class IpLineSearchBase
    {
    public:
        /**
         * @brief Find an acceptable trial point.
         */
        virtual void find_acceptable_trial_point() = 0;

        /**
         * @brief Reset the line search algorithm.
         */
        virtual void reset() = 0;

        /**
         * @brief Reset the line search for a new iteration.
         */
        virtual void reset_linesearch() = 0;

        virtual void register_options(OptionRegistry &registry) = 0;
        virtual bool is_acceptable_to_current_iterate(const Scalar trial_barr,
                                                      const Scalar trial_theta,
                                                      const bool called_from_resto = false) = 0;
        virtual bool is_acceptable_to_filter(const Scalar trial_barr, const Scalar trial_theta) = 0;

    protected:
        virtual ~IpLineSearchBase() = default;
    };
    /**
     *
     * todo: do we really want to template on the linear solver type here?
     *
     * */

    /**
     * @brief Concrete implementation of line search for a specific problem type.
     *
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename LinearSolverType, typename ProblemType>
    class IpLinesearch : public IpLineSearchBase
    {
        typedef std::shared_ptr<LinearSolverType> LinearSolverSp;
        typedef std::shared_ptr<IpData<ProblemType>> IpDataSp;
        typedef std::shared_ptr<Nlp<ProblemType>> NlpSp;
        typedef IpIterate<ProblemType> IpIterateType;
        typedef std::shared_ptr<IpRestoPhaseBase> RestorationPhaseSp;

    public:
        /**
         * @brief Construct a new IpLinesearch object.
         *
         * @param ipdata Shared pointer to the interior point algorithm data.
         * @param linear_solver Shared pointer to the linear solver.
         * @param restoration_phase Shared pointer to the restoration phase.
         */
        IpLinesearch(const IpDataSp &ipdata, const LinearSolverSp &linear_solver,
                     const RestorationPhaseSp &restoration_phase);

        void find_acceptable_trial_point() override;
        void reset() override;
        void reset_linesearch() override;
        // Setter methods for options
        void set_max_soc(const Index &value) { max_soc_ = value; }
        void set_kappa_soc(const Scalar &value) { kappa_soc_ = value; }
        void set_theta_min(const Scalar &value) { theta_min_ = value; }
        void set_theta_max(const Scalar &value) { theta_max_ = value; }
        void set_theta_min_fact(const Scalar &value) { theta_min_fact_ = value; }
        void set_s_phi(const Scalar &value) { s_phi_ = value; }
        void set_s_theta(const Scalar &value) { s_theta_ = value; }
        void set_delta(const Scalar &value) { delta_ = value; }
        void set_tiny_step_tol(const Scalar &value) { tiny_step_tol_ = value; }
        void set_tiny_step_y_tol(const Scalar &value) { tiny_step_y_tol_ = value; }
        void set_gamma_theta(const Scalar &value) { gamma_theta_ = value; }
        void set_gamma_phi(const Scalar &value) { gamma_phi_ = value; }
        void set_eta_phi(const Scalar &value) { eta_phi_ = value; }
        void set_alpha_min_frac(const Scalar &value) { alpha_min_frac_ = value; }
        void set_theta_max_fact(const Scalar &value) { theta_max_fact_ = value; }
        void set_max_filter_resets(const Index &value) { max_filter_resets_ = value; }
        void set_filter_reset_trigger(const Index &value) { filter_reset_trigger_ = value; }
        void set_obj_max_incr(const Scalar &value) { obj_max_incr_ = value; }
        void set_watchdog_shortened_iter_trigger(const Index &value)
        {
            watchdog_shortened_iter_trigger_ = value;
        }
        void set_watchdog_trial_iter_max(const Index &value) { watchdog_trial_iter_max_ = value; }
        void set_alpha_red_factor(const Scalar &value) { alpha_red_factor_ = value; }
        void set_max_iter(const Index &max_iter)
        {
            filter().reserve(max_iter_ + 1);
            max_iter_ = max_iter;
        }
        void set_soft_rest_pd_error_reduction_factor(const Scalar &value)
        {
            soft_rest_pd_error_reduction_factor_ = value;
        }
        void set_max_soft_resto_iters(const Index &value) { max_soft_resto_iters_ = value; }
        void register_options(OptionRegistry &registry) override;

    private:
        IpFilter &filter();
        void init_this_line_search(bool in_watchdog);
        bool do_backtracking_line_search(bool skip_first_trial_point, Scalar &alpha_primal,
                                         bool &corr_taken, bool &soc_taken, Index &n_steps,
                                         bool &evaluation_error);
        void start_watchdog();
        void stop_watchdog();
        bool check_acceptability_trial_point(const Scalar alpha_primal);
        void perform_dual_step(const Scalar alpha_primal, const Scalar alpha_dual);
        bool try_second_order_correction(const Scalar alpha_primal_test, Scalar &alpha_primal);
        bool try_soft_resto_step(bool &satisfies_original_criterion);
        bool is_acceptable_to_current_iterate(const Scalar trial_barr, const Scalar trial_theta,
                                              const bool called_from_resto = false) override;
        bool is_acceptable_to_filter(const Scalar trial_barr, const Scalar trial_theta) override;
        bool armijo_holds(const Scalar alpha_primal);
        bool detect_tiny_step();
        bool is_f_type(const Scalar alpha_primal);
        bool check_acceptability_of_trial_point(const Scalar alpha_primal);
        void augment_filter();
        Scalar compute_alpha_min();
        char update_for_next_iteration(const Scalar alpha_primal_test);
        void update_step_info(const Scalar alpha_primal, const Scalar alpha_dual,
                              const Index n_steps, const char info_alpha_primal_char);
        void prepare_resto_phase_start();

        IpFilter filter_;
        IpDataSp ipdata_;
        LinearSolverSp linear_solver_;
        RestorationPhaseSp restoration_phase_;
        VecRealAllocated soc_rhs_x_;
        VecRealAllocated soc_rhs_s_;
        VecRealAllocated soc_rhs_g_;
        VecRealAllocated soc_g_accumulate_;
        VecRealAllocated soc_rhs_cl_;
        VecRealAllocated soc_rhs_cu_;
        // options
        Index max_soc_ = 4;
        Scalar kappa_soc_ = 0.99;
        Scalar theta_min_ = -1.;
        Scalar theta_max_ = -1.;
        Scalar theta_min_fact_ = 1e-4;
        Scalar s_phi_ = 2.3;
        Scalar s_theta_ = 1.1;
        Scalar delta_ = 1.0;
        Scalar tiny_step_tol_ = 10 * std::numeric_limits<Scalar>::epsilon();
        Scalar tiny_step_y_tol_ = 1e-2;
        Scalar gamma_theta_ = 1e-5;
        Scalar gamma_phi_ = 1e-8;
        Scalar eta_phi_ = 1e-8;
        Scalar alpha_min_frac_ = 0.05;
        Scalar theta_max_fact_ = 1e4;
        Index max_filter_resets_ = 5;
        Index filter_reset_trigger_ = 5;
        Scalar obj_max_incr_ = 5.;
        Index watchdog_shortened_iter_trigger_ = 10;
        Index watchdog_trial_iter_max_ = 3;
        Scalar alpha_red_factor_ = 0.5;
        Index max_iter_ = 1000;
        Scalar soft_rest_pd_error_reduction_factor_ = 1. - 1e-4;
        Index max_soft_resto_iters_ = 10;

        // internal staticstics
        bool in_watchdog_;
        Scalar reference_theta_;
        Scalar reference_barr_;
        Scalar reference_grad_bar_delta_;
        Scalar watchdog_theta_;
        Scalar watchdog_barr_;
        Scalar watchdog_grad_bar_delta_;
        Scalar watchdog_alpha_primal_test_;
        bool last_rejection_due_to_filter_;
        Index filter_reset_count_;
        Index filter_reject_count_;
        bool tiny_step_last_iteration_;
        Index sucessive_tiny_step_count_;
        Index acceptable_iteration_count_;
        Scalar last_mu_;
        Index watchdog_shortened_iter_count_;
        Index watchdog_trial_iter_count_ = 0;
        bool in_soft_resto_phase_;
        Index soft_resto_counter_;
    };
} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_linesearch_hpp__
