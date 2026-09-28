// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_IP_HESSIAN_EVALUATOR_HPP
#define FATROP_IP_HESSIAN_EVALUATOR_HPP

#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/nlp/hessian.hpp"
#include "fatrop/common/options.hpp"
#include <memory>

namespace fatrop
{
    /// Hessian allocation, evaluation and history lifecycle, independent of the
    /// iterate's derivative cache. The default evaluator uses the NLP callback.
    template <typename P> class IpHessianEvaluator
    {
    public:
        virtual ~IpHessianEvaluator() = default;
        virtual std::unique_ptr<Hessian<P>> create(const ProblemDims<P> &dims) const
        {
            return std::make_unique<Hessian<P>>(dims);
        }
        virtual void evaluate(IpIterate<P> &point, IpIterate<P> &current, Hessian<P> &out)
        {
            const Index status = point.nlp()->eval_lag_hess(
                point.info(), point.objective_scale, point.primal_x(), point.primal_s(),
                point.dual_eq(), out);
            fatrop_assert_msg(status == 0, "Error in evaluating the Hessian of the Lagrangian.");
        }
        virtual bool uses_history() const { return false; }
        virtual void reset() {}
        virtual void register_options(OptionRegistry &) {}
    };
}
#endif
