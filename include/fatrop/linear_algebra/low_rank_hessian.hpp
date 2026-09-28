// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_LOW_RANK_HESSIAN_HPP
#define FATROP_LOW_RANK_HESSIAN_HPP

#include "fatrop/linear_algebra/matrix.hpp"
#include <memory>

namespace fatrop
{
    /** B = sigma I + W diag(+1,-1,+1,-1,...) W^T.
     * The structured Hessian stores sigma I; only the correction lives here.
     * Immutable snapshots let trial and watchdog iterates retain their Hessian.
     */
    struct LowRankHessian
    {
        LowRankHessian(Index n, Index pairs, Scalar scale)
            : sigma(scale), W(n, 2 * pairs), pairs_(pairs) { W = 0.; }

        Scalar sigma;
        MatRealAllocated W;
        Index rank() const { return 2 * pairs_; }
        void reset(Index pairs, Scalar scale)
        {
            fatrop_assert_msg(2 * pairs <= W.n(), "Low-rank Hessian capacity exceeded");
            pairs_ = pairs;
            sigma = scale;
        }
        static Scalar sign(Index j) { return j % 2 == 0 ? 1. : -1.; }

        void add_product(const VecRealView &x, VecRealView &out) const
        {
            for (Index j = 0; j < rank(); ++j)
            {
                Scalar a = 0.;
                for (Index i = 0; i < W.m(); ++i) a += W(i, j) * x(i);
                a *= sign(j);
                for (Index i = 0; i < W.m(); ++i) out(i) += a * W(i, j);
            }
        }
    private:
        Index pairs_;
    };
}
#endif
