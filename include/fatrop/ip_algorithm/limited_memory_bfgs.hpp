// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_LIMITED_MEMORY_BFGS_HPP
#define FATROP_LIMITED_MEMORY_BFGS_HPP

#include "fatrop/common/options.hpp"
#include "fatrop/linear_algebra/low_rank_hessian.hpp"
#include "fatrop/nlp/jacobian.hpp"
#include <algorithm>
#include <cmath>
#include <array>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace fatrop
{
    /** Standard Euclidean L-BFGS of the Lagrangian, shared by all NLP backends.
     * Secant differences use the NEW multipliers at BOTH endpoints. Historical
     * pairs are frozen; only one previous Jacobian and objective gradient is kept.
     */
    template <typename ProblemType> class LimitedMemoryBfgs
    {
    public:
        explicit LimitedMemoryBfgs(const ProblemInfo<ProblemType> &info, Index capacity = 6)
            : capacity_(checked_capacity(capacity)), max_history_(std::min(Index(6), capacity_)),
              S_(info.number_of_tangent_variables, capacity_),
              Y_(info.number_of_tangent_variables, capacity_),
              previous_x_(info.number_of_tangent_variables),
              previous_gradient_(info.number_of_tangent_variables),
              old_lag_(info.number_of_tangent_variables), new_lag_(info.number_of_tangent_variables),
              step_(info.number_of_tangent_variables), difference_(info.number_of_tangent_variables),
              bs_(info.number_of_tangent_variables), previous_jacobian_(info.dims)
        {
            for (auto &slot : snapshots_)
                slot = std::make_shared<LowRankHessian>(info.number_of_tangent_variables, capacity_, 1.);
        }

        bool enabled() const { return enabled_; }
        Index history_size() const { return count_; }
        Index skipped_updates() const { return skipped_; }

        void reset()
        {
            count_ = start_ = skipped_ = 0;
            previous_valid_ = false;
            current_.reset();
            sigma_ = init_val_;
        }

        std::shared_ptr<const LowRankHessian> update(
            const ProblemInfo<ProblemType> &info, const VecRealView &x,
            const VecRealView &gradient, const Jacobian<ProblemType> &jacobian,
            const VecRealView &lambda)
        {
            const Index n = gradient.m();
            if (x.m() != n || n != S_.m())
                throw std::invalid_argument("limited-memory requires matching Euclidean primal coordinates");
            if (!previous_valid_ || skipped_ >= max_skipping_)
            {
                count_ = start_ = skipped_ = 0;
                sigma_ = init_val_;
                snapshot(info, x, gradient, jacobian);
                current_ = acquire_snapshot(0, sigma_);
                return current_;
            }

            old_lag_ = previous_gradient_;
            previous_jacobian_.transpose_apply_on_right(info, lambda, 1., old_lag_, old_lag_);
            new_lag_ = gradient;
            jacobian.transpose_apply_on_right(info, lambda, 1., new_lag_, new_lag_);
            Scalar ss = 0., sy = 0., yy = 0., smax = 0.;
            for (Index i = 0; i < n; ++i)
            {
                step_(i) = x(i) - previous_x_(i);
                difference_(i) = new_lag_(i) - old_lag_(i);
                ss += step_(i) * step_(i);
                sy += step_(i) * difference_(i);
                yy += difference_(i) * difference_(i);
                smax = std::max(smax, std::abs(step_(i)));
            }
            const Scalar eps = std::numeric_limits<Scalar>::epsilon();
            if (smax < 100. * eps || !std::isfinite(ss) || !std::isfinite(sy) ||
                !std::isfinite(yy) || sy <= std::sqrt(eps) * std::sqrt(ss) * std::sqrt(yy))
            {
                ++skipped_;
                snapshot(info, x, gradient, jacobian);
                return current_;
            }

            sigma_ = initialization_ == Initialization::Constant ? init_val_ :
                     initialization_ == Initialization::Scalar2 ? yy / sy : sy / ss;
            sigma_ = std::clamp(sigma_, init_min_, init_max_);
            if (count_ == max_history_) start_ = (start_ + 1) % capacity_;
            else ++count_;
            const Index newest = (start_ + count_ - 1) % capacity_;
            for (Index i = 0; i < n; ++i)
            {
                S_(i, newest) = step_(i);
                Y_(i, newest) = difference_(i);
            }

            // Replay retained BFGS pairs into an unused preallocated snapshot.
            // The current, trial, stored and factored models keep their own snapshots.
            auto next = acquire_snapshot(count_, sigma_);
            bool valid = true;
            for (Index j = 0; j < count_; ++j)
            {
                const Index pair = (start_ + j) % capacity_;
                for (Index i = 0; i < n; ++i) bs_(i) = sigma_ * S_(i, pair);
                for (Index k = 0; k < 2 * j; ++k)
                {
                    Scalar a = 0.;
                    for (Index i = 0; i < n; ++i) a += next->W(i, k) * S_(i, pair);
                    a *= LowRankHessian::sign(k);
                    for (Index i = 0; i < n; ++i) bs_(i) += a * next->W(i, k);
                }
                Scalar sbs = 0., sty = 0.;
                for (Index i = 0; i < n; ++i)
                {
                    sbs += S_(i, pair) * bs_(i);
                    sty += S_(i, pair) * Y_(i, pair);
                }
                if (!(sbs > 0.) || !std::isfinite(sbs)) { valid = false; break; }
                for (Index i = 0; i < n; ++i)
                {
                    next->W(i, 2 * j) = Y_(i, pair) / std::sqrt(sty);
                    next->W(i, 2 * j + 1) = bs_(i) / std::sqrt(sbs);
                }
            }
            if (!valid)
            {
                count_ = start_ = 0;
                sigma_ = init_val_;
                next->reset(0, sigma_);
            }
            current_ = std::move(next);
            skipped_ = 0;
            snapshot(info, x, gradient, jacobian);
            return current_;
        }

        void set_approximation(const std::string &value)
        {
            if (value != "exact" && value != "limited-memory")
                throw std::invalid_argument("hessian_approximation must be exact or limited-memory");
            enabled_ = value == "limited-memory";
            reset();
        }
        void set_history(const Index &value)
        {
            if (value < 1 || value > capacity_)
                throw std::invalid_argument("limited_memory_max_history exceeds construction capacity or is nonpositive");
            max_history_ = value;
            reset();
        }
        void set_max_skipping(const Index &value)
        {
            if (value < 1) throw std::invalid_argument("limited_memory_max_skipping must be positive");
            max_skipping_ = value;
        }
        void set_initialization(const std::string &value)
        {
            if (value != "scalar1" && value != "scalar2" && value != "constant")
                throw std::invalid_argument("limited_memory_initialization must be scalar1, scalar2 or constant");
            initialization_ = value == "scalar1" ? Initialization::Scalar1 :
                              value == "scalar2" ? Initialization::Scalar2 : Initialization::Constant;
        }
        void set_init_val(const Scalar &value) { positive(value); init_val_ = value; }
        void set_init_min(const Scalar &value)
        {
            positive(value);
            if (value > init_max_) throw std::invalid_argument("limited_memory_init_val_min exceeds maximum");
            init_min_ = value;
        }
        void set_init_max(const Scalar &value)
        {
            positive(value);
            if (value < init_min_) throw std::invalid_argument("limited_memory_init_val_max is below minimum");
            init_max_ = value;
        }
        void register_options(OptionRegistry &registry)
        {
            registry.register_option("hessian_approximation", &LimitedMemoryBfgs::set_approximation, this);
            registry.register_option("limited_memory_max_history", &LimitedMemoryBfgs::set_history, this);
            registry.register_option("limited_memory_max_skipping", &LimitedMemoryBfgs::set_max_skipping, this);
            registry.register_option("limited_memory_initialization", &LimitedMemoryBfgs::set_initialization, this);
            registry.register_option("limited_memory_init_val", &LimitedMemoryBfgs::set_init_val, this);
            registry.register_option("limited_memory_init_val_min", &LimitedMemoryBfgs::set_init_min, this);
            registry.register_option("limited_memory_init_val_max", &LimitedMemoryBfgs::set_init_max, this);
        }

    private:
        enum class Initialization { Scalar1, Scalar2, Constant };
        const Index capacity_;
        bool enabled_ = false, previous_valid_ = false;
        Index max_history_, max_skipping_ = 2, skipped_ = 0, count_ = 0, start_ = 0;
        Scalar init_val_ = 1., init_min_ = 1e-8, init_max_ = 1e8, sigma_ = 1.;
        Initialization initialization_ = Initialization::Scalar1;
        MatRealAllocated S_, Y_;
        VecRealAllocated previous_x_, previous_gradient_, old_lag_, new_lag_, step_, difference_, bs_;
        Jacobian<ProblemType> previous_jacobian_;
        // At most five distinct snapshots are held by the history, three iterates,
        // and the factored solver. The sixth is the next writable snapshot.
        std::array<std::shared_ptr<LowRankHessian>, 6> snapshots_;
        std::shared_ptr<const LowRankHessian> current_;

        static Index checked_capacity(Index capacity)
        {
            if (capacity < 1) throw std::invalid_argument("Hessian history capacity must be positive");
            return capacity;
        }
        std::shared_ptr<LowRankHessian> acquire_snapshot(Index pairs, Scalar scale)
        {
            for (auto &slot : snapshots_)
                if (slot.use_count() == 1)
                {
                    slot->reset(pairs, scale);
                    return slot;
                }
            throw std::logic_error("Too many retained limited-memory Hessian snapshots");
        }

        static void positive(Scalar value)
        {
            if (!(value > 0.) || !std::isfinite(value))
                throw std::invalid_argument("limited-memory scale must be finite and positive");
        }
        void snapshot(const ProblemInfo<ProblemType> &info, const VecRealView &x,
                      const VecRealView &gradient, const Jacobian<ProblemType> &jacobian)
        {
            previous_jacobian_.copy_from(jacobian);
            previous_x_ = x;
            previous_gradient_ = gradient;
            previous_valid_ = true;
        }
    };
}
#endif
