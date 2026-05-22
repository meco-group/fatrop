//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_dense_jacobian_hpp__
#define __fatrop_dense_jacobian_hpp__

#include "fatrop/dense/fwd.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/nlp/jacobian.hpp"
#include <iosfwd>

namespace fatrop
{
    typedef ProblemInfo<DenseType> DenseInfo;

    /**
     * @brief Specialization of the constraint Jacobian for dense problems.
     *
     * Stored in transposed form (one extra row reserved for the right-hand side),
     * matching the convention used for the OCP variant.
     */
    template <> struct Jacobian<DenseType>
    {
        Jacobian(const ProblemDims<DenseType> &dims);

        /**
         * @brief Equality-constraint Jacobian (transposed).
         *
         * Dimensions: (nx + 1) x ng.
         *  - Rows [0:nx) hold the Jacobian columns of g_eq w.r.t. x.
         *  - Row nx holds the right-hand side g_eq(x_k) (last row).
         */
        MatRealAllocated Gg_eqt;

        /**
         * @brief Inequality-constraint Jacobian (transposed).
         *
         * Dimensions: (nx + 1) x ng_ineq. Same layout as @c Gg_eqt.
         */
        MatRealAllocated Gg_ineqt;

        void apply_on_right(const DenseInfo &info, const VecRealView &x, Scalar alpha,
                            const VecRealView &y, VecRealView &out) const;
        void transpose_apply_on_right(const DenseInfo &info, const VecRealView &mult_eq,
                                      Scalar alpha, const VecRealView &y, VecRealView &out) const;
        void get_rhs(const DenseInfo &info, VecRealView &rhs) const;
        void set_rhs(const DenseInfo &info, const VecRealView &rhs);

        friend std::ostream &operator<<(std::ostream &os, const Jacobian &jac);
    };
} // namespace fatrop

#endif // __fatrop_dense_jacobian_hpp__
