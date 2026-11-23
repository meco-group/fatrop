//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_graph_blockwise_abstract_hpp__
#define __fatrop_graph_graph_blockwise_abstract_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/graph/graph_abstract.hpp"
#include "fatrop/linear_algebra/fwd.hpp"

#include <memory>
#include <utility>
#include <vector>

namespace fatrop
{
    /**
     * @brief Block-by-block interface for a graph-structured NLP.
     *
     * Same problem as @ref GraphAbstract, evaluated one block at a time;
     * wrap it in a @ref GraphBlockwiseAdapter to solve it.
     *
     * The Hessian / gradient and Jacobian evaluation hooks are invoked once
     * per block per iteration: @c eval_Hk fills the (i, j) lower-triangle
     * block of the Hessian, zeroed on entry; @c eval_grad_k
     * fills block @c k of the gradient; @c eval_Ggt_ineq_k fills the
     * @c (n_k + 1) x ng_ineq[k] transposed Jacobian + RHS row for block @c k.
     */
    class GraphBlockwiseAbstract
    {
    public:
        /// Number of variable blocks @c N.
        virtual Index get_num_blocks() const = 0;

        /// Size of the primal-variable block @c k (manifold dimension).
        virtual Index get_block_size(Index k) const = 0;

        /// Tangent-space dimension of block @c k. Defaults to @c get_block_size(k).
        /// Override only when block @c k lives on a Lie group / manifold and
        /// the search direction needs fewer components than the primal —
        /// e.g. a unit quaternion with @c get_block_size(k) == 4 and
        /// @c get_block_size_tan(k) == 3 (so(3)).
        virtual Index get_block_size_tan(Index k) const { return get_block_size(k); }

        /// Number of block-local inequality constraints for block @c k.
        /// Defaults to 0 (no inequality constraints).
        virtual Index get_ng_ineq(Index) const { return 0; }

        /// Off-diagonal Hessian edges (unordered pairs, duplicates allowed).
        virtual std::vector<std::pair<Index, Index>> get_off_diag_edges() const = 0;

        // ---------------------------------------------------------------
        // Per-iteration evaluations.
        // ---------------------------------------------------------------

        /**
         * @brief Compute one Hessian block.
         *
         * Writes the (@c i, @c j) lower-triangle Hessian block (@c i >= @c j)
         * of the Lagrangian, with the convention that @c x[k] points to the
         * primal variables of block @c k and @c mult[k] points to the
         * inequality-constraint multipliers of block @c k.
         *
         * @param i, j           Block indices in the lower triangle (@c i >= @c j).
         * @param objective_scale Scalar multiplier on the objective.
         * @param x              Flat primal vector (all blocks concatenated).
         * @param mult           Flat inequality-multiplier vector
         *                       (length sum of @c get_ng_ineq).
         * @param res            BLASFEO matrix view for the block to fill,
         *                       shape @c get_block_size(i) x @c get_block_size(j).
         */
        virtual Index eval_Hk(Index i, Index j, const Scalar *objective_scale, const Scalar *x,
                              const Scalar *mult, MAT *res) = 0;

        /**
         * @brief Compute the @c k-th block of the objective gradient.
         *
         * @param k              Block index.
         * @param objective_scale Scalar multiplier on the objective.
         * @param x              Flat primal vector (all blocks concatenated).
         * @param res            Output vector of size @c get_block_size(k).
         */
        virtual Index eval_grad_k(Index k, const Scalar *objective_scale, const Scalar *x,
                                  Scalar *res) = 0;

        /**
         * @brief Compute the transposed inequality Jacobian + RHS row for
         *        block @c k.
         *
         * The matrix @c res has shape @c (n_k + 1) x ng_ineq[k]. The top
         * @c n_k rows hold @c d g_ineq_k / d x_k (in tangent space). The
         * trailing row (index @c n_k) stores @c g_ineq_k(x_k).
         */
        virtual Index eval_Ggt_ineq_k(Index, const Scalar *, MAT *) { return 0; }

        /**
         * @brief Compute the inequality constraint values @c g_ineq_k(x_k)
         *        into a flat output vector positioned at the per-block
         *        offset.
         *
         * Default: no-op (zero inequality constraints).
         */
        virtual Index eval_gineq_k(Index, const Scalar *, Scalar *) { return 0; }

        /// Evaluate the objective.
        virtual Index eval_f(const Scalar *objective_scale, const Scalar *x, Scalar *res) = 0;

        /**
         * @brief Lower / upper bounds for the inequality constraints.
         *
         * Layout: per-block concatenation. Use @c +inf / @c -inf for
         * unbounded sides. Default: no-op (zero inequality constraints).
         */
        virtual Index get_bounds(Scalar *, Scalar *) const { return 0; }

        /// Initial primal guess (flat vector, all blocks concatenated).
        virtual Index get_initial(Scalar *x) const = 0;

        /**
         * @brief Per-block retraction hook for manifold-valued primals.
         *
         * Computes @c x_next[k] = retract_k(x[k], alpha * delta_x[k]) for each
         * block @c k. The default implementation is the Euclidean update
         * @c x_next = x + alpha * delta_x — valid whenever every block is
         * Euclidean (@c get_block_size_tan(k) == get_block_size(k) for all
         * @c k). Override to support Lie-group blocks (e.g. SE(2), SO(3),
         * unit quaternions).
         *
         * @param x        Flat primal vector (sum of @c get_block_size).
         * @param delta_x  Flat tangent vector (sum of @c get_block_size_tan).
         * @param alpha    Step size.
         * @param x_next   Output flat primal vector (sum of @c get_block_size).
         */
        virtual void apply_retraction(const Scalar *x, const Scalar *delta_x, const Scalar alpha,
                                      Scalar *x_next)
        {
            // Euclidean fallback. Requires get_block_size_tan == get_block_size
            // for every block.
            Index off_primal = 0;
            Index off_tan = 0;
            const Index N = get_num_blocks();
            for (Index k = 0; k < N; ++k)
            {
                const Index n = get_block_size(k);
                for (Index i = 0; i < n; ++i)
                    x_next[off_primal + i] = x[off_primal + i] + alpha * delta_x[off_tan + i];
                off_primal += n;
                off_tan += get_block_size_tan(k);
            }
        }

        virtual ~GraphBlockwiseAbstract() = default;
    };

    /**
     * @brief @ref GraphAbstract that evaluates a @ref GraphBlockwiseAbstract
     *        block by block.
     */
    class GraphBlockwiseAdapter : public GraphAbstract
    {
    public:
        GraphBlockwiseAdapter(const std::shared_ptr<GraphBlockwiseAbstract> &graph);

        GraphStructure get_structure() const override { return structure_; }
        Index eval_hessian(const Scalar *objective_scale, const Scalar *x, const Scalar *mult,
                           GraphHessianView &H) override;
        Index eval_gradient(const Scalar *objective_scale, const Scalar *x,
                            Scalar *grad) override;
        Index eval_jacobian_ineq(const Scalar *x, GraphJacobianView &J) override;
        Index eval_gineq(const Scalar *x, Scalar *g) override;
        Index eval_f(const Scalar *objective_scale, const Scalar *x, Scalar *res) override
        {
            return graph_->eval_f(objective_scale, x, res);
        }
        Index get_bounds(Scalar *lower, Scalar *upper) const override
        {
            return graph_->get_bounds(lower, upper);
        }
        Index get_initial(Scalar *x) const override { return graph_->get_initial(x); }
        void apply_retraction(const Scalar *x, const Scalar *delta_x, const Scalar alpha,
                              Scalar *x_next) override
        {
            graph_->apply_retraction(x, delta_x, alpha, x_next);
        }

    private:
        std::shared_ptr<GraphBlockwiseAbstract> graph_;
        GraphStructure structure_; // off_diag_edges normalised to unique (i > j)
        std::vector<Index> offsets_tan_;
        std::vector<Index> offsets_ineq_;
    };
} // namespace fatrop

#endif // __fatrop_graph_graph_blockwise_abstract_hpp__
