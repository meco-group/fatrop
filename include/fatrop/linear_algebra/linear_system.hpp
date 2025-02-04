//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

/**
 * @file linear_system.hpp
 * @brief Defines the LinearSystem class template for representing and manipulating linear systems.
 *
 * This header file provides the LinearSystem class template, which represents a linear system
 * of the form Ax = -b. It offers methods for accessing and modifying the right-hand side (RHS)
 * of the equation and applying the system matrix to a vector.
 */

#ifndef __fatrop_linear_algebra_linear_system_hpp__
#define __fatrop_linear_algebra_linear_system_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/fwd.hpp"

namespace fatrop
{
    /**
     * @brief Represents a linear system of the form Ax = -b.
     *
     * @tparam LsType The underlying type used to represent the linear system.
     */
    template <typename LsType> class LinearSystem;
    // {
    // public:
    //     /**
    //      * @brief Get the number of rows in the linear system.
    //      *
    //      * @return Index The number of rows.
    //      */
    //     Index m();

    //     /**
    //      * @brief Get the right-hand side (RHS) of the linear system.
    //      *
    //      * @param[out] out VecRealView to store the RHS.
    //      */
    //     void get_rhs(VecRealView &out);

    //     /**
    //      * @brief Set the right-hand side (RHS) of the linear system.
    //      *
    //      * @param[in] in VecRealView containing the new RHS values.
    //      */
    //     void set_rhs(const VecRealView &in);

    //     /**
    //      * @brief Apply the system matrix A to a vector x on the right (i.e., compute Ax).
    //      *
    //      * @param[in] x VecRealView representing the input vector.
    //      * @param[out] out VecRealView to store the result of Ax.
    //      */
    //     void apply_on_right(const VecRealView &x, Scalar alpha, const VecRealView& y, VecRealView
    //     &out);
    // };
} // namespace fatrop

#endif // __fatrop_linear_algebra_linear_system_hpp__
