//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_linesearch_hpp__
#define __fatrop_ip_algorithm_ip_linesearch_hpp__
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

        /**
         * @brief Accept the current trial iterate as the new iterate.
         */
        virtual void accept_trial_iterate() = 0;

    protected:
        virtual ~IpLineSearchBase() = default;
    };

    /**
     * @brief Concrete implementation of line search for a specific problem type.
     *
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename ProblemType> class IpLinesearch : public IpLineSearchBase
    {
        typedef std::shared_ptr<PdSolverOrig<ProblemType>> PdSolverSp;
        typedef std::shared_ptr<IpData<ProblemType>> IpDataSp;
        typedef std::shared_ptr<Nlp<ProblemType>> NlpSp;
        typedef IpIterate<ProblemType> IpIterateType;

    public:
        /**
         * @brief Construct a new IpLinesearch object.
         *
         * @param ipdata Shared pointer to the interior point algorithm data.
         * @param linear_solver Shared pointer to the primal-dual linear solver.
         */
        IpLinesearch(const IpDataSp &ipdata, const PdSolverSp &linear_solver);

        void find_acceptable_trial_point() override;
        void reset() override;
        void reset_linesearch() override;
        void accept_trial_iterate() override;

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
        bool is_acceptable_to_current_iterate(const Scalar trial_barr, const Scalar trial_theta,
                                              const bool called_from_resto = false);
        bool is_acceptable_to_filter(const Scalar trial_barr, const Scalar trial_theta);
        bool armijo_holds(const Scalar alpha_primal);
        bool detect_tiny_step();
        bool is_f_type(const Scalar alpha_primal);
        bool check_acceptability_of_trial_point(const Scalar alpha_primal);
        void augment_filter();
        Scalar compute_alpha_min();
        char update_for_next_iteration(const Scalar alpha_primal_test);
        void update_step_info(const Scalar alpha_primal, const Scalar alpha_dual, const Index n_steps, const char info_alpha_primal_char);

        IpFilter filter_;
        IpDataSp ipdata_;
        PdSolverSp linear_solver_;
        VecRealAllocated soc_rhs_x_;
        VecRealAllocated soc_rhs_s_;
        VecRealAllocated soc_rhs_g_;
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
    };
} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_linesearch_hpp__
