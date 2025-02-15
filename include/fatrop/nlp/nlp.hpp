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
    // forward declaration
    template <typename ProblemType> class IpData;
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
                                              const VecRealView &primal_x,
                                              const VecRealView &primal_s, VecRealView &grad_x,
                                              VecRealView &grad_s) = 0;
        virtual Index eval_objective(const ProblemInfo<ProblemType> &info,
                                     const Scalar objective_scale, const VecRealView &primal_x,
                                     const VecRealView &primal_s, Scalar &res) = 0;
        virtual Index get_bounds(const ProblemInfo<ProblemType> &info, VecRealView &lower_bounds,
                                 VecRealView &upper_bounds) = 0;
        virtual Index get_initial_primal(const ProblemInfo<ProblemType> &info,
                                         VecRealView &primal_x) = 0;
        virtual void get_primal_damping(const ProblemInfo<ProblemType> &info,
                                        VecRealView &damping) = 0;
        virtual void apply_jacobian_s_transpose(const ProblemInfo<ProblemType> &info,
                                                const VecRealView &multipliers, const Scalar alpha,
                                                const VecRealView &y, VecRealView &out) = 0;
        virtual void callback(const IpData<ProblemType> &ip_data) {};
        virtual ~Nlp() = default;
    };
} // namespace fatrop

#endif //__fatrop_nlp_nlp__
