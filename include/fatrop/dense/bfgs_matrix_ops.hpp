// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_DENSE_BFGS_MATRIX_OPS_HPP
#define FATROP_DENSE_BFGS_MATRIX_OPS_HPP

#include "fatrop/dense/hessian.hpp"
#include "fatrop/nlp/bfgs_matrix_ops.hpp"

namespace fatrop
{
    template <> struct BfgsMatrixOps<DenseType>
    {
        template <typename F> static void for_each_block(const Hessian<DenseType> &h, F &&f)
        { f(0, 0, h.Hht.n()); }

        template <typename F> static void for_each_lower(
            const Hessian<DenseType> &h, Index, Index, F &&f)
        {
            for (Index i = 0; i < h.Hht.n(); ++i)
                for (Index j = 0; j <= i; ++j) f(i, j, h.Hht(i, j));
        }
        static void set(Hessian<DenseType> &h, Index, Index, Index i, Index j, Scalar value)
        {
            h.Hht(i, j) = value;
            h.Hht(j, i) = value;
        }
        static void stabilize(Hessian<DenseType> &, VecRealView &) {}
        static void copy(const Hessian<DenseType> &source, Hessian<DenseType> &target)
        { target.Hht = source.Hht; }
    };
}
#endif
