// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_LIMITED_MEMORY_HESSIAN_EVALUATOR_HPP
#define FATROP_LIMITED_MEMORY_HESSIAN_EVALUATOR_HPP

#include "fatrop/ip_algorithm/ip_hessian_evaluator.hpp"
#include "fatrop/ip_algorithm/limited_memory_bfgs.hpp"
#include "fatrop/nlp/hessian_with_low_rank.hpp"

namespace fatrop
{
    template <typename P> class LimitedMemoryHessianEvaluator : public IpHessianEvaluator<P>
    {
    public:
        explicit LimitedMemoryHessianEvaluator(const ProblemInfo<P> &info, Index history_capacity = 6)
            : bfgs_(info, history_capacity) {}
        std::unique_ptr<Hessian<P>> create(const ProblemDims<P> &dims) const override
        {
            return std::make_unique<HessianWithLowRank<P>>(dims);
        }
        void evaluate(IpIterate<P> &point, IpIterate<P> &current, Hessian<P> &out) override
        {
            auto &extended = static_cast<HessianWithLowRank<P> &>(out);
            if (!bfgs_.enabled())
            {
                extended.low_rank.reset();
                IpHessianEvaluator<P>::evaluate(point, current, out);
                return;
            }
            // Trials reuse the accepted model. Only an accepted point forms a
            // secant pair, always using its multipliers at both endpoints.
            auto approximation = &point == &current
                ? bfgs_.update(point.info(), point.primal_x(), point.obj_gradient_x(),
                               point.jacobian(), point.dual_eq())
                : static_cast<HessianWithLowRank<P> &>(current.hessian()).low_rank;
            out.set_diagonal(approximation->sigma);
            extended.low_rank = std::move(approximation);
        }
        bool uses_history() const override { return bfgs_.enabled(); }
        void reset() override { bfgs_.reset(); }
        void register_options(OptionRegistry &registry) override { bfgs_.register_options(registry); }
        // The selector owns hessian_approximation, while this evaluator owns
        // the limited-memory-specific settings.
        void register_memory_options(OptionRegistry &registry) { bfgs_.register_options(registry, false); }
        void set_approximation(const std::string &value) { bfgs_.set_approximation(value); }
        const LimitedMemoryBfgs<P> &approximation() const { return bfgs_; }

    private:
        LimitedMemoryBfgs<P> bfgs_;
    };
}
#endif
