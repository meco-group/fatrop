//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//
/**
 * @file ip_filter.hpp
 * @brief Defines the IpFilter class for interior-point methods.
 */

#ifndef __fatrop_ip_algorithm_ip_filter_hpp__
#define __fatrop_ip_algorithm_ip_filter_hpp__
#include "fatrop/context/context.hpp"
#include <vector>

namespace fatrop
{
    /**
     * @brief Data structure for filter entries.
     */
    struct IpFilterData
    {
        /**
         * @brief Objective function value.
         */
        Scalar obj_value;
        /**
         * @brief Constraint violation.
         */
        Scalar constr_viol;
    };

    /**
     * @brief IpFilter class for interior-point methods.
     *
     * This class implements an filter for interior-point methods, as described in Ipopt
     * implementation paper.
     */
    class IpFilter
    {
    public:
        /**
         * @brief Constructor.
         */
        IpFilter();
        /**
         * @brief Reset the filter.
         */
        void reset();
        /**
         * @brief Reserve space for the filter.
         * @param[in] max_size Maximum size of the filter.
         */
        void reserve(const Index max_size);
        /**
         * @brief Check if a data point is acceptable to the filter.
         * @param[in] data Data point to check.
         * @return True if the data point is acceptable, false otherwise.
         */
        bool is_acceptable(const IpFilterData &data) const;
        /**
         * @brief Augment the filter with a new data point.
         * @param[in] data Data point to augment the filter with.
         */
        void augment(const IpFilterData &data);
        /**
         * @brief Get the current size of the filter.
         * @return Current size of the filter.
         */
        Index size() const;

    private:
        /**
         * @brief Check if data point a dominates data point b.
         * @param[in] a Data point a.
         * @param[in] b Data point b.
         * @return True if a dominates b, false otherwise.
         */
        static bool a_dominates_b(const IpFilterData &a, const IpFilterData &b);
        /**
         * @brief Vector of filter data points.
         */
        std::vector<IpFilterData> data_;
    };
} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_algorithm_hpp__
