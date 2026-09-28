// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_FULL_BFGS_HESSIAN_EVALUATOR_HPP
#define FATROP_FULL_BFGS_HESSIAN_EVALUATOR_HPP

#include "fatrop/ip_algorithm/full_bfgs.hpp"
#include "fatrop/ip_algorithm/ip_hessian_evaluator.hpp"
#include "fatrop/nlp/hessian_with_low_rank.hpp"

namespace fatrop
{
    template <typename P> class FullBfgsHessianEvaluator : public IpHessianEvaluator<P>
    {
    public:
        explicit FullBfgsHessianEvaluator(const ProblemInfo<P> &info) : bfgs_(info) {}

        std::unique_ptr<Hessian<P>> create(const ProblemDims<P> &dims) const override
        {
            return std::make_unique<HessianWithLowRank<P>>(dims);
        }

        void evaluate(IpIterate<P> &point, IpIterate<P> &current, Hessian<P> &out) override
        {
            auto &extended = static_cast<HessianWithLowRank<P> &>(out);
            extended.low_rank.reset();
            if (&point == &current)
                bfgs_.update(point.info(), point.primal_x(), point.obj_gradient_x(),
                             point.jacobian(), point.dual_eq());
            bfgs_.copy_to(out);
        }

        bool uses_history() const override { return true; }
        void reset() override { bfgs_.reset(); }

    private:
        FullBfgs<P> bfgs_;
    };
}
#endif
