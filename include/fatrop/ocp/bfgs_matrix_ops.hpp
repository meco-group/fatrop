// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_OCP_BFGS_MATRIX_OPS_HPP
#define FATROP_OCP_BFGS_MATRIX_OPS_HPP

#include "fatrop/ocp/hessian.hpp"
#include "fatrop/nlp/bfgs_matrix_ops.hpp"

namespace fatrop
{
    template <> struct BfgsMatrixOps<OcpType>
    {
        template <typename F> static void for_each_block(const Hessian<OcpType> &h, F &&f)
        {
            Index offset = 0;
            for (Index k = 0; k < static_cast<Index>(h.RSQrqt.size()); ++k)
            {
                f(k, offset, h.RSQrqt[k].n());
                offset += h.RSQrqt[k].n();
            }
        }
        template <typename F> static void for_each_lower(
            const Hessian<OcpType> &h, Index k, Index offset, F &&f)
        {
            const auto &stage = h.RSQrqt[k];
            for (Index i = 0; i < stage.n(); ++i)
                for (Index j = 0; j <= i; ++j)
                    f(offset + i, offset + j, stage(i, j));
        }
        static void set(Hessian<OcpType> &h, Index k, Index offset,
                        Index i, Index j, Scalar value)
        {
            auto &stage = h.RSQrqt[k];
            stage(i - offset, j - offset) = value;
            stage(j - offset, i - offset) = value;
        }
        static void stabilize(Hessian<OcpType> &, VecRealView &) {}
        static void copy(const Hessian<OcpType> &source, Hessian<OcpType> &target)
        {
            for (Index k = 0; k < static_cast<Index>(source.RSQrqt.size()); ++k)
                target.RSQrqt[k] = source.RSQrqt[k];
        }
    };
}
#endif
