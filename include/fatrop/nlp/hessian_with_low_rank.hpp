// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_HESSIAN_WITH_LOW_RANK_HPP
#define FATROP_HESSIAN_WITH_LOW_RANK_HPP

#include "fatrop/nlp/hessian.hpp"
#include "fatrop/linear_algebra/low_rank_hessian.hpp"

namespace fatrop
{
    /// Optional operator extension; native Hessians retain their original storage.
    template <typename P> class HessianWithLowRank : public Hessian<P>
    {
    public:
        using Hessian<P>::Hessian;
        std::shared_ptr<const LowRankHessian> low_rank;

        void apply_on_right(const ProblemInfo<P> &info, const VecRealView &x,
                            Scalar alpha, const VecRealView &y, VecRealView &out) const override
        {
            Hessian<P>::apply_on_right(info, x, alpha, y, out);
            if (low_rank) low_rank->add_product(x, out);
        }

        void set_zero() override
        {
            low_rank.reset();
            Hessian<P>::set_zero();
        }
    };
}
#endif
