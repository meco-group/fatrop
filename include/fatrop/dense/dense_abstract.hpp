//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_dense_dense_abstract_hpp__
#define __fatrop_dense_dense_abstract_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/fwd.hpp"

namespace fatrop
{
    /**
     * @brief User-facing abstract interface for a dense NLP problem.
     *
     * Dense problem:
     *     min   f(x)
     *     s.t.  g_eq(x)    = 0
     *           L <= g_ineq(x) <= U
     *
     * The primal variable @c x can either be a plain Euclidean vector or live
     * on a manifold (Lie group). In the manifold case the search direction
     * @c delta_x lives in the tangent space and is recovered through the
     * @c apply_retraction hook. The size of @c x is @c get_nx(); the size of
     * @c delta_x is @c get_nx_tangent() (defaults to @c get_nx() for the
     * Euclidean case). The Hessian and Jacobian buffers and the gradient row
     * passed by fatrop are sized in tangent space:
     *
     *  - @c eval_Hh writes a (nx_tan + 1) x nx_tan matrix whose top nx_tan block
     *    is the Hessian of the Lagrangian (in tangent space) and whose last row
     *    is the gradient.
     *  - @c eval_Ggt writes a (nx_tan + 1) x ng matrix; top nx_tan rows are
     *    @c grad_x(g_eq)^T (tangent-space), last row is @c g_eq(x).
     *  - @c eval_Ggt_ineq is the same for inequality constraints.
     *
     * The lam vector passed to @c eval_Hh has length @c ng + @c ng_ineq with
     * the equality multipliers followed by the slack-equality (inequality)
     * multipliers.
     */
    class DenseAbstract
    {
    public:
        virtual Index get_nx() const = 0;
        virtual Index get_ng() const = 0;
        virtual Index get_ng_ineq() const = 0;

        /**
         * @brief Tangent-space dimension of the primal variable.
         *
         * Defaults to @c get_nx(). Override only when @c x lives on a manifold
         * (Lie-group optimisation) and the search direction needs fewer
         * components than the primal — for example a unit quaternion with
         * @c get_nx() == 4 and @c get_nx_tangent() == 3 (so(3)).
         */
        virtual Index get_nx_tangent() const { return get_nx(); }

        virtual Index eval_Hh(const Scalar *objective_scale, const Scalar *x, const Scalar *lam,
                              MAT *res) = 0;
        virtual Index eval_Ggt(const Scalar *x, MAT *res) = 0;
        virtual Index eval_Ggt_ineq(const Scalar *x, MAT *res) = 0;
        virtual Index eval_g(const Scalar *x, Scalar *res) = 0;
        virtual Index eval_gineq(const Scalar *x, Scalar *res) = 0;
        virtual Index eval_grad(const Scalar *objective_scale, const Scalar *x, Scalar *res) = 0;
        virtual Index eval_f(const Scalar *objective_scale, const Scalar *x, Scalar *res) = 0;
        virtual Index get_bounds(Scalar *lower, Scalar *upper) const = 0;
        virtual Index get_initial(Scalar *x) const = 0;

        /**
         * @brief Retraction hook for manifold-valued primals.
         *
         * Computes @c x_next = retract(x, alpha * delta_x). The default is the
         * Euclidean update @c x_next = x + alpha * delta_x (works whenever
         * @c get_nx_tangent() == get_nx()).
         *
         * @param x       current primal, size @c get_nx()
         * @param delta_x search direction, size @c get_nx_tangent()
         * @param alpha   step size
         * @param x_next  retracted primal, size @c get_nx() (output)
         */
        virtual void apply_retraction(const Scalar *x, const Scalar *delta_x, const Scalar alpha,
                                      Scalar *x_next)
        {
            const Index n = get_nx();
            for (Index i = 0; i < n; ++i)
                x_next[i] = x[i] + alpha * delta_x[i];
        }

        /**
         * @brief Map the dual of the *scaled* equality constraints to the
         *        dual of the unscaled constraints before @c eval_Hh sees it.
         *
         * If @c eval_Ggt or @c eval_g pre-scaled the linearised constraints by
         * some matrix @c M(x), fatrop's Newton solver returns
         * @c lambda_tilde of the scaled constraint and this hook recovers
         * @c lambda = M(x)^T lambda_tilde. The default is the identity.
         *
         * @param x         current primal, size @c get_nx()
         * @param dual_in   dual of the scaled constraint, length ng + ng_ineq
         * @param dual_out  dual of the unscaled constraint, length ng + ng_ineq
         */
        virtual void apply_dual_eq_transformation(const Scalar *x, const Scalar *dual_in,
                                                  Scalar *dual_out)
        {
            const Index n = get_ng() + get_ng_ineq();
            for (Index i = 0; i < n; ++i)
                dual_out[i] = dual_in[i];
        }

        virtual ~DenseAbstract() = default;
    };
} // namespace fatrop

#endif // __fatrop_dense_dense_abstract_hpp__
