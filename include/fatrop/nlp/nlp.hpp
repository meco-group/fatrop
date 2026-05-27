//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_nlp_nlp__
#define __fatrop_nlp_nlp__

#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"
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
        /**
         * @brief Compute the retracted primal variable
         *        primal_x_next = retract(primal_x, alpha * delta_primal_x).
         *
         * The default implementation is the Euclidean update
         * @code primal_x_next = primal_x + alpha * delta_primal_x; @endcode
         * which requires @c number_of_variables == @c number_of_tangent_variables.
         *
         * Override this method to support Lie-group / manifold optimization where
         * @c primal_x lives on a manifold and @c delta_primal_x lives in the corresponding
         * tangent / Lie algebra. The retraction must satisfy @c retract(x, 0) == x .
         *
         * Dimensions:
         *  - @c primal_x and @c primal_x_next have size @c nlp_dims().number_of_variables
         *  - @c delta_primal_x has size @c nlp_dims().number_of_tangent_variables
         */
        virtual void apply_retraction(const ProblemInfo<ProblemType> &info,
                                      const VecRealView &primal_x,
                                      const VecRealView &delta_primal_x, const Scalar alpha,
                                      VecRealView &primal_x_next)
        {
            primal_x_next = primal_x + alpha * delta_primal_x;
        }
        virtual void callback(const IpData<ProblemType> &ip_data) {};
        virtual ~Nlp() = default;
    };
} // namespace fatrop

#endif //__fatrop_nlp_nlp__
