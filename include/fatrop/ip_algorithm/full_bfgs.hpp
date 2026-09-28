// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_FULL_BFGS_HPP
#define FATROP_FULL_BFGS_HPP

#include "fatrop/nlp/bfgs_matrix_ops.hpp"
#include "fatrop/nlp/jacobian.hpp"
#include "fatrop/nlp/hessian.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fatrop
{
    /** Explicit BFGS secant history and update formula. Matrix storage and
     * sparsity are supplied by the backend's BfgsMatrixOps specialization.
     */
    template <typename P> class FullBfgs
    {
    public:
        explicit FullBfgs(const ProblemInfo<P> &info)
            : model_(info.dims), previous_x_(info.number_of_tangent_variables),
              previous_gradient_(info.number_of_tangent_variables),
              old_lag_(info.number_of_tangent_variables), new_lag_(info.number_of_tangent_variables),
              step_(info.number_of_tangent_variables), difference_(info.number_of_tangent_variables),
              bs_(info.number_of_tangent_variables), row_sum_(info.number_of_tangent_variables),
              previous_jacobian_(info.dims)
        { reset(); }

        void reset() { previous_valid_ = false; model_.set_diagonal(1.); }

        void update(const ProblemInfo<P> &info, const VecRealView &x,
                    const VecRealView &gradient, const Jacobian<P> &jacobian,
                    const VecRealView &lambda)
        {
            const Index n = gradient.m();
            if (x.m() != n || n != step_.m())
                throw std::invalid_argument("BFGS requires matching Euclidean primal coordinates");
            if (previous_valid_)
            {
                old_lag_ = previous_gradient_;
                previous_jacobian_.transpose_apply_on_right(info, lambda, 1., old_lag_, old_lag_);
                new_lag_ = gradient;
                jacobian.transpose_apply_on_right(info, lambda, 1., new_lag_, new_lag_);
                for (Index i = 0; i < n; ++i)
                {
                    step_(i) = x(i) - previous_x_(i);
                    difference_(i) = new_lag_(i) - old_lag_(i);
                }
                BfgsMatrixOps<P>::for_each_block(model_, [&](Index block, Index offset, Index length) {
                    update_block(block, offset, length);
                });
            }
            previous_x_ = x;
            previous_gradient_ = gradient;
            previous_jacobian_.copy_from(jacobian);
            previous_valid_ = true;
        }

        void copy_to(Hessian<P> &out) const { BfgsMatrixOps<P>::copy(model_, out); }

    private:
        Hessian<P> model_;
        VecRealAllocated previous_x_, previous_gradient_, old_lag_, new_lag_;
        VecRealAllocated step_, difference_, bs_, row_sum_;
        Jacobian<P> previous_jacobian_;
        bool previous_valid_ = false;

        void update_block(Index block, Index offset, Index length)
        {
            Scalar ss = 0., sy = 0., yy = 0., smax = 0.;
            for (Index i = offset; i < offset + length; ++i)
            {
                ss += step_(i) * step_(i);
                sy += step_(i) * difference_(i);
                yy += difference_(i) * difference_(i);
                smax = std::max(smax, std::abs(step_(i)));
                bs_(i) = 0.;
            }
            const Scalar eps = std::numeric_limits<Scalar>::epsilon();
            if (smax < 100. * eps || !std::isfinite(ss) || !std::isfinite(sy) ||
                !std::isfinite(yy) || sy <= std::sqrt(eps) * std::sqrt(ss * yy)) return;

            BfgsMatrixOps<P>::for_each_lower(model_, block, offset, [&](Index i, Index j, Scalar v) {
                bs_(i) += v * step_(j);
                if (i != j) bs_(j) += v * step_(i);
            });
            Scalar sbs = 0.;
            for (Index i = offset; i < offset + length; ++i) sbs += step_(i) * bs_(i);
            if (!(sbs > 0.) || !std::isfinite(sbs)) { model_.set_diagonal(1.); return; }

            BfgsMatrixOps<P>::for_each_lower(model_, block, offset, [&](Index i, Index j, Scalar v) {
                const Scalar updated = v + difference_(i) * difference_(j) / sy
                                           - bs_(i) * bs_(j) / sbs;
                BfgsMatrixOps<P>::set(model_, block, offset, i, j, updated);
            });
            BfgsMatrixOps<P>::stabilize(model_, row_sum_);
        }
    };
}
#endif
