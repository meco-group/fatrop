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
     * as well as other important algorithm parameters like the barrier parameter (mu)
     * and the iteration count.
     *
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename ProblemType> struct IpData
    {
        typedef std::shared_ptr<Nlp<ProblemType>> NlpSp;
        typedef IpIterate<ProblemType> Iterate;

        /**
         * @brief Construct a new IpData object.
         *
         * @param nlp Shared pointer to the NLP problem.
         */
        IpData(const NlpSp &nlp);

        void reset();

        /**
         * @brief Accept the trial iterate as the new current iterate.
         *
         * This method switches the trial iterate to become the current iterate
         * and resets the trial iterate.
         */
        void accept_trial_iterate();

        /**
         * @brief Backup the current iterate into the stored iterate.
         *
         * This method is used for the watchdog mechanism. It works by switching
         * the current and stored iterate pointers. Be careful: the current iterate
         * is invalidated after this operation and becomes meaningless.
         */
        void backup_current_iterate();

        /**
         * @brief Restore the stored iterate as the current iterate.
         *
         * This method is used in conjunction with backup_current_iterate()
         * to implement the watchdog mechanism.
         */
        void restore_current_iterate();

        /**
         * @brief Get a reference to the current iterate.
         *
         * @return Iterate& Reference to the current iterate.
         */
        Iterate &current_iterate();

        /**
         * @brief Get a const reference to the current iterate.
         *
         * @return const Iterate& Const reference to the current iterate.
         */
        const Iterate &current_iterate() const;
        /**
         * @brief Get a reference to the trial iterate.
         *
         * @return Iterate& Reference to the trial iterate.
         */
        Iterate &trial_iterate() { return *trial_iterate_; }

        /**
         * @brief Get a const reference to the trial iterate.
         *
         * @return const Iterate& Const reference to the trial iterate.
         */
        const Iterate &trial_iterate() const { return *trial_iterate_; }

        /**
         * @brief Get a reference to the stored iterate.
         *
         * @return Iterate& Reference to the stored iterate.
         */
        Iterate &stored_iterate() { return *stored_iterate_; }

        /**
         * @brief Get a const reference to the stored iterate.
         *
         * @return const Iterate& Const reference to the stored iterate.
         */
        const Iterate &stored_iterate() const { return *stored_iterate_; }

        /**
         * @brief Set the barrier parameter mu.
         *
         * @param mu The new value for mu.
         */
        void set_mu(const Scalar mu);

        /**
         * @brief Get the current value of the barrier parameter mu.
         *
         * @return Scalar The current value of mu.
         */
        Scalar mu() const { return mu_; }

        /**
         * @brief Get the current iteration number.
         *
         * @return Index The current iteration number.
         */
        Index iteration_number() const { return iteration_number_; }

        /**
         * @brief Set the iteration number.
         *
         * @param iteration_number The new iteration number.
         */
        void set_iteration_number(const Index iteration_number)
        {
            iteration_number_ = iteration_number;
        }

        /**
         * @brief Invalidate the current iterate.
         *
         * This method marks the current iterate as invalid.
         */
        void invalidate_current_iterate() { current_iterate_is_valid_ = false; }

        /**
         * @brief Validate the current iterate.
         *
         * This method marks the current iterate as valid.
         */
        void validate_current_iterate() { current_iterate_is_valid_ = true; }

    private:
        Scalar mu_;                ///< Barrier value of the NLP.
        Index iteration_number_;   ///< Number of the current iteration.
        Iterate iterate_data_[3];  ///< Data for the three iterates (current, trial, and stored).
        Iterate *current_iterate_; ///< Pointer to the current iterate.
        Iterate *trial_iterate_;   ///< Pointer to the trial iterate.
        Iterate *stored_iterate_;  ///< Pointer to the stored iterate.
        bool current_iterate_is_valid_ = true; ///< Flag indicating if the current iterate is valid.
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
