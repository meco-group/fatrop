// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_LIMITED_MEMORY_NLP_RESTO_HPP
#define FATROP_LIMITED_MEMORY_NLP_RESTO_HPP

#include "fatrop/ip_algorithm/ip_nlp_resto.hpp"
#include "fatrop/ip_algorithm/ip_hessian_evaluator.hpp"
#include "fatrop/ip_algorithm/limited_memory_hessian_evaluator.hpp"

namespace fatrop
{
    /// Approximate Hessians include the restoration objective's curvature in x. Its slack
    /// curvature remains explicit because the approximation only acts on x.
    template <typename P> class LimitedMemoryNlpResto : public IpNlpResto<P>
    {
    public:
        LimitedMemoryNlpResto(
            const std::shared_ptr<Nlp<P>> &nlp,
            std::shared_ptr<IpHessianEvaluator<P>> evaluator)
            : IpNlpResto<P>(nlp), evaluator_(std::move(evaluator)) {}

        void get_primal_damping(const ProblemInfo<P> &info, VecRealView &damping) override
        {
            if (!evaluator_->uses_history())
            {
                IpNlpResto<P>::get_primal_damping(info, damping);
                return;
            }
            auto xs = damping.block(info.number_of_primal_variables +
                                    info.number_of_slack_variables, 0);
            this->nlp_orig_->get_primal_damping(info, xs);
            auto ds = xs.block(info.number_of_slack_variables, info.offset_slack);
            auto drs = this->dr_.block(info.number_of_slack_variables, info.offset_slack);
            ds = ds + this->zeta_ * drs * drs;
            damping.block(info.number_of_eq_constraints, info.offset_slack_p) = 0.;
            damping.block(info.number_of_eq_constraints, info.offset_slack_n) = 0.;
        }

    private:
        std::shared_ptr<IpHessianEvaluator<P>> evaluator_;
    };
}
#endif
