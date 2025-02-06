//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_data_hpp__
#define __fatrop_ip_algorithm_ip_data_hpp__
#include "fatrop/common/exception.hpp"
#include "fatrop/context/context.hpp"
#include "fatrop/ip_algorithm/ip_iterate.hpp"
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

        /**
         * @brief Construct a new IpData object.
         * @param nlp Shared pointer to the NLP problem.
         */
        IpData(const NlpSp &nlp);

        /**
         * @brief Reset the IpData object to its initial state.
         */
        void reset();

        /**
         * @brief Get the NLP problem associated with this data.
         * @return NlpSp The shared pointer to the NLP problem.
         */
        NlpSp get_nlp() const;

        /**
         * @brief Accept the trial iterate as the new current iterate.
         * Switches the trial iterate to become the current iterate and resets the trial iterate.
         */
        void accept_trial_iterate();

        /**
         * @brief Backup the current iterate into the stored iterate.
         * Used for the watchdog mechanism. Switches current and stored iterate pointers.
         * Warning: Invalidates the current iterate.
         */
        void backup_current_iterate();

        /**
         * @brief Restore the stored iterate as the current iterate.
         * Used in conjunction with backup_current_iterate() for the watchdog mechanism.
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
         * @brief Invalidate the current iterate.
         * Marks the current iterate as invalid.
         */
        void invalidate_current_iterate() { current_iterate_is_valid_ = false; }

        /**
         * @brief Validate the current iterate.
         * Marks the current iterate as valid.
         */
        void validate_current_iterate() { current_iterate_is_valid_ = true; }

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

    private:
        Index iteration_number_;   ///< Number of the current iteration.
        Iterate iterate_data_[3];  ///< Data for the three iterates (current, trial, and stored).
        Iterate *current_iterate_; ///< Pointer to the current iterate.
        Iterate *trial_iterate_;   ///< Pointer to the trial iterate.
        Iterate *stored_iterate_;  ///< Pointer to the stored iterate.
        bool current_iterate_is_valid_ = true; ///< Flag indicating if the current iterate is valid.
        bool tiny_step_flag_ = false;
        Scalar tol_ = 1e-8;
    };
    template <typename ProblemType> IpIterate<ProblemType> &IpData<ProblemType>::current_iterate()
    {
        fatrop_assert_msg(current_iterate_is_valid_, "the current iterate is invalidated.");
        return *current_iterate_;
    }
    template <typename ProblemType>
    const IpIterate<ProblemType> &IpData<ProblemType>::current_iterate() const
    {
        fatrop_assert_msg(current_iterate_is_valid_, "the current iterate is invalidated.");
        return *current_iterate_;
    }

} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_data_hpp__
