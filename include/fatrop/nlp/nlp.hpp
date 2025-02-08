//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_nlp_nlp__
#define __fatrop_nlp_nlp__

#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/nlp/fwd.hpp"

namespace fatrop
{
    template <typename ProblemType> class Nlp
    {
    public:
        virtual const NlpDims &nlp_dims() const = 0;
        virtual const ProblemDims<ProblemType> &problem_dims() const = 0;
        virtual Index eval_lag_hess(const ProblemInfo<ProblemType> &info,
                                    const Scalar objective_scale, const VecRealView &primal_x,
                                    const VecRealView &primal_s, const VecRealView &lam,
                                    Hessian<ProblemType> &hess) = 0;
        virtual Index eval_constr_jac(const ProblemInfo<ProblemType> &info,
                                      const VecRealView &primal_x, const VecRealView &primal_s,
                                      Jacobian<ProblemType> &jac) = 0;
        virtual Index eval_constraint_violation(const ProblemInfo<ProblemType> &info,
                                                const VecRealView &primal_x,
                                                const VecRealView &primal_s, VecRealView &res) = 0;
        virtual Index eval_objective_gradient(const ProblemInfo<ProblemType> &info,
                                              const Scalar objective_scale,
                                              const VecRealView &primal_x, VecRealView &grad_x,
                                              VecRealView &grad_s) = 0;
        virtual Index eval_objective(const ProblemInfo<ProblemType> &info,
                                     const Scalar objective_scale, const VecRealView &primal_x,
                                     const VecRealView &primal_s, Scalar &res) = 0;
        virtual Index get_bounds(const ProblemInfo<ProblemType> &info, VecRealView &lower_bounds,
                                 VecRealView &upper_bounds) = 0;
        virtual Index get_initial_primal(const ProblemInfo<ProblemType> &info, VecRealView &primal_x) = 0;
        virtual ~Nlp() = default;
    };
} // namespace fatrop

#endif //__fatrop_nlp_nlp__
