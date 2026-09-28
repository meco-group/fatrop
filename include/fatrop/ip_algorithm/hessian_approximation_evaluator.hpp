// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_HESSIAN_APPROXIMATION_EVALUATOR_HPP
#define FATROP_HESSIAN_APPROXIMATION_EVALUATOR_HPP

#include "fatrop/ip_algorithm/full_bfgs_hessian_evaluator.hpp"
#include "fatrop/ip_algorithm/limited_memory_hessian_evaluator.hpp"

namespace fatrop
{
    /** Selects the Hessian evaluator without putting approximation state in the
     * algorithm, native Hessians, or linear solvers. The limited-memory evaluator
     * also handles the exact callback when its approximation is disabled.
     */
    template <typename P> class HessianApproximationEvaluator : public IpHessianEvaluator<P>
    {
    public:
        HessianApproximationEvaluator(const ProblemInfo<P> &info, Index history_capacity = 6)
            : limited_(info, history_capacity), full_(info), selected_(&limited_) {}

        std::unique_ptr<Hessian<P>> create(const ProblemDims<P> &dims) const override
        {
            return limited_.create(dims);
        }
        void evaluate(IpIterate<P> &point, IpIterate<P> &current, Hessian<P> &out) override
        {
            selected_->evaluate(point, current, out);
        }
        bool uses_history() const override { return selected_->uses_history(); }
        void reset() override
        {
            limited_.reset();
            full_.reset();
        }
        void register_options(OptionRegistry &registry) override
        {
            registry.register_option("hessian_approximation", &HessianApproximationEvaluator::set_approximation, this);
            limited_.register_memory_options(registry);
        }
        void set_approximation(const std::string &value)
        {
            if (value != "exact" && value != "limited-memory" && value != "bfgs")
                throw std::invalid_argument("hessian_approximation must be exact, limited-memory or bfgs");
            limited_.set_approximation(value == "limited-memory" ? value : "exact");
            full_.reset();
            selected_ = value == "bfgs" ? static_cast<IpHessianEvaluator<P> *>(&full_)
                                        : static_cast<IpHessianEvaluator<P> *>(&limited_);
        }

        const LimitedMemoryBfgs<P> &limited_memory() const { return limited_.approximation(); }

    private:
        LimitedMemoryHessianEvaluator<P> limited_;
        FullBfgsHessianEvaluator<P> full_;
        IpHessianEvaluator<P> *selected_;
    };
}
#endif
