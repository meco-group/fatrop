//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_data_hpp__
#define __fatrop_ip_algorithm_ip_data_hpp__
#include "fatrop/common/exception.hpp"
#include "fatrop/common/fwd.hpp"
#include "fatrop/context/context.hpp"
#include "fatrop/ip_algorithm/ip_iterate.hpp"
#include "fatrop/ip_algorithm/ip_timings.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include <vector>

namespace fatrop
{
    /**
     * @brief Stores and manages data for the interior point algorithm.
     *
     * This class maintains the current, trial, and stored iterates of the algorithm,
     * as well as other important algorithm parameters like the iteration count and tolerance.
     *
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename ProblemType> struct IpData
    {
        typedef std::shared_ptr<Nlp<ProblemType>> NlpSp;
        typedef IpIterate<ProblemType> Iterate;
        typedef ProblemInfo<ProblemType> InfoType;

        /**
         * @brief Construct a new IpData object.
         * @param nlp Shared pointer to the NLP problem.
         */
        IpData(const NlpSp &nlp);

        const InfoType &info() { return info_; }

        /**
         * @brief Reset the IpData object to its initial state.
         */
        void reset();

        /**
         * @brief Get the NLP problem associated with this data.
         * @return NlpSp The shared pointer to the NLP problem.
         */
        NlpSp get_nlp() const { return nlp_; };

        /**
         * @brief Accept the trial iterate as the new current iterate.
         * Switches the trial iterate to become the current iterate and resets the trial iterate.
         */
        void accept_trial_iterate();

        /**
         * @brief Backup the current iterate into the stored iterate.
         * Used for the watchdog mechanism. Switches current and stored iterate pointers.
         */
        void store_current_iterate();

        /**
         * @brief Restore the stored iterate as the current iterate.
         * Used in conjunction with store_current_iterate() for the watchdog mechanism.
         */
        void restore_current_iterate();

        /**
         * @brief Get a reference to the current iterate.
         * @return Iterate& Reference to the current iterate.
         */
        Iterate &current_iterate();

        /**
         * @brief Get a const reference to the current iterate.
         * @return const Iterate& Const reference to the current iterate.
         */
        const Iterate &current_iterate() const;

        /**
         * @brief Get a reference to the trial iterate.
         * @return Iterate& Reference to the trial iterate.
         */
        Iterate &trial_iterate() { return *trial_iterate_; }

        /**
         * @brief Get a const reference to the trial iterate.
         * @return const Iterate& Const reference to the trial iterate.
         */
        const Iterate &trial_iterate() const { return *trial_iterate_; }

        /**
         * @brief Get a reference to the stored iterate.
         * @return Iterate& Reference to the stored iterate.
         */
        Iterate &stored_iterate() { return *stored_iterate_; }

        /**
         * @brief Get a const reference to the stored iterate.
         * @return const Iterate& Const reference to the stored iterate.
         */
        const Iterate &stored_iterate() const { return *stored_iterate_; }

        /**
         * @brief Get the current iteration number.
         * @return Index The current iteration number.
         */
        Index iteration_number() const { return iteration_number_; }

        /**
         * @brief Set the iteration number.
         * @param iteration_number The new iteration number.
         */
        void set_iteration_number(const Index iteration_number)
        {
            iteration_number_ = iteration_number;
        }

        /**
         * @brief Check if a tiny step flag is set.
         * @return bool True if the tiny step flag is set, false otherwise.
         */
        bool tiny_step_flag() const { return tiny_step_flag_; }

        /**
         * @brief Set the tiny step flag.
         * @param tiny_step_flag The new value for the tiny step flag.
         */
        void set_tiny_step_flag(bool tiny_step_flag) { tiny_step_flag_ = tiny_step_flag; }

        /**
         * @brief Get the current tolerance value.
         * @return Scalar The current tolerance value.
         */
        Scalar tolerance() const { return tol_; }

        // Setter method for tolerance
        void set_tolerance(const Scalar &value) { tol_ = value; }

        // Register options
        void register_options(OptionRegistry &registry);

        IpTimingStatistics &timing_statistics() { return timings_; }

    private:
        NlpSp nlp_; ///< Shared pointer to the NLP problem.
        InfoType info_;
        Index iteration_number_;   ///< Number of the current iteration.
        Iterate iterate_data_[3];  ///< Data for the three iterates (current, trial, and stored).
        Iterate *current_iterate_; ///< Pointer to the current iterate.
        Iterate *trial_iterate_;   ///< Pointer to the trial iterate.
        Iterate *stored_iterate_;  ///< Pointer to the stored iterate.
        bool tiny_step_flag_ = false;
        Scalar tol_ = 1e-8;
        IpTimingStatistics timings_;
        Hessian<ProblemType> hessian_data_[2];
        Jacobian<ProblemType> jacobian_data_[2];
        Hessian<ProblemType> *hessian_curr_ = nullptr;
        Jacobian<ProblemType> *jacobian_curr_ = nullptr;
        Hessian<ProblemType> *hessian_stored_ = nullptr;
        Jacobian<ProblemType> *jacobian_stored_ = nullptr;
        bool stored_iterate_is_valid_ = false;
        // problem information
    protected:
        friend class IpIterate<ProblemType>;
        VecRealAllocated lower_bounds_; ///< Lower bounds of the variables.
        VecRealAllocated upper_bounds_; ///< Upper bounds of the variables.
        std::vector<bool>
            lower_bounded_; ///< Boolean vector indicating if the variables are lower bounded.
        std::vector<bool>
            upper_bounded_; ///< Boolean vector indicating if the variables are upper bounded.
        std::vector<bool> single_lower_bounded_; ///< Boolean vector indicating if the variables are
                                                 ///< lower bounded.
        std::vector<bool> single_upper_bounded_; ///< Boolean vector indicating if the variables are
                                                 ///< upper bounded.
        Index number_of_bounds_ ; ///< Total number of bounds in the problem
    };
    template <typename ProblemType> IpIterate<ProblemType> &IpData<ProblemType>::current_iterate()
    {
        return *current_iterate_;
    }
    template <typename ProblemType>
    const IpIterate<ProblemType> &IpData<ProblemType>::current_iterate() const
    {
        return *current_iterate_;
    }

} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_data_hpp__
