//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_dense_hessian_hpp__
#define __fatrop_dense_hessian_hpp__

#include "fatrop/dense/fwd.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/nlp/hessian.hpp"
#include <iosfwd>

namespace fatrop
{
    typedef ProblemInfo<DenseType> DenseInfo;

    /**
     * @brief Specialization of the Lagrangian Hessian for dense problems.
     *
     * Stored as a single (nx + 1) x nx block. The top nx-by-nx block is the
     * Hessian @c H; the last row holds the objective gradient @c h (so the
     * matrix is @c [H; h^T]^T — hence the name @c Hht). This is the layout the
     * augmented-system solver and the BLASFEO kernels expect.
     */
    template <> struct Hessian<DenseType>
    {
        Hessian(const ProblemDims<DenseType> &dims);

        /// Hessian + gradient block, shape (nx + 1) x nx; last row is the gradient.
        MatRealAllocated Hht;

        void apply_on_right(const DenseInfo &info, const VecRealView &x, Scalar alpha,
                            const VecRealView &y, VecRealView &out) const;
        void get_rhs(const DenseInfo &info, VecRealView &out) const;
        void set_rhs(const DenseInfo &info, const VecRealView &in);
        void set_zero();

        friend std::ostream &operator<<(std::ostream &os, const Hessian<DenseType> &hess);
    };
} // namespace fatrop

#endif // __fatrop_dense_hessian_hpp__
