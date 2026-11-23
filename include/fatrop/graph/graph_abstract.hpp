//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_graph_abstract_hpp__
#define __fatrop_graph_graph_abstract_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/graph/block_pd_matrix.hpp"
#include "fatrop/linear_algebra/matrix.hpp"

#include <utility>
#include <vector>

namespace fatrop
{
    /**
     * @brief Structure of a graph-structured NLP.
     *
     * @c block_sizes_tan defaults to @c block_sizes and @c ng_ineq to zeros when
     * left empty. @c off_diag_edges lists every structurally nonzero
     * off-diagonal Hessian block once, as @c (i, j) with @c i > j.
     */
    struct GraphStructure
    {
        std::vector<Index> block_sizes;
        std::vector<Index> block_sizes_tan;
        std::vector<Index> ng_ineq;
        std::vector<std::pair<Index, Index>> off_diag_edges;
    };

    /**
     * @brief Writable view on the Lagrangian Hessian and objective gradient.
     *
     *   - @c diag(k): block @c (k, k), @c n_k_tan x n_k_tan; only the lower
     *     triangle is read.
     *   - @c off(e): block @c (i, j) for @c (i, j) = off_diag_edges[e],
     *     @c n_i_tan x n_j_tan (lower triangle, @c i > j).
     *   - @c grad(): flat tangent-space gradient (block_offsets_tan layout).
     *
     * All blocks and the gradient are zero on entry.
     */
    class GraphHessianView
    {
    public:
        GraphHessianView(BlockPdMatrix &H, const Index *diag_ids, const Index *off_ids,
                         Scalar *grad)
            : H_(H), diag_ids_(diag_ids), off_ids_(off_ids), grad_(grad)
        {
        }
        MatRealView diag(Index k) { return H_.block_by_id(diag_ids_[k]); }
        MatRealView off(Index e) { return H_.block_by_id(off_ids_[e]); }
        Scalar *grad() { return grad_; }

    private:
        BlockPdMatrix &H_;
        const Index *diag_ids_;
        const Index *off_ids_;
        Scalar *grad_;
    };

    /**
     * @brief Writable view on the per-block inequality Jacobians.
     *
     * @c block(k) is the @c (n_k_tan + 1) x ng_ineq[k] transposed Jacobian of
     * @c g_ineq_k w.r.t. @c x_k (tangent space), with @c g_ineq_k(x_k) in the
     * trailing row. All blocks are zero on entry.
     */
    class GraphJacobianView
    {
    public:
        GraphJacobianView(MatRealAllocated *blocks) : blocks_(blocks) {}
        MatRealView block(Index k) { return blocks_[k]; }

    private:
        MatRealAllocated *blocks_;
    };

    /**
     * @brief User-facing abstract interface for a graph-structured NLP.
     *
     * Graph problem:
     *     min   f(x_0, ..., x_{N-1})
     *     s.t.  L_k_i <= g_k_i(x_k) <= U_k_i    for each block k, constraint i
     *
     * where the primal variable @c x is partitioned into @c N blocks
     * @c x_0, ..., x_{N-1}. The Lagrangian Hessian has a symmetric
     * block-sparse structure described by @c GraphStructure::off_diag_edges;
     * inequality constraints are block-local (each depends on exactly one
     * block of primal variables); equality constraints are not supported.
     *
     * Index conventions follow the rest of fatrop: @c x is the flat
     * concatenation of the primal blocks, gradients and search directions the
     * flat concatenation of the tangent blocks, and inequality vectors the
     * concatenation of the per-block inequality vectors.
     *
     * Every evaluation covers the whole problem in one call.
     */
    class GraphAbstract
    {
    public:
        /// Block sizes, inequality counts and Hessian sparsity. Queried once.
        virtual GraphStructure get_structure() const = 0;

        /**
         * @brief Lagrangian Hessian and objective gradient.
         *
         * @param objective_scale Scalar multiplier on the objective.
         * @param x               Flat primal vector.
         * @param mult            Flat inequality-multiplier vector.
         * @param H               Hessian blocks and gradient to fill.
         */
        virtual Index eval_hessian(const Scalar *objective_scale, const Scalar *x,
                                   const Scalar *mult, GraphHessianView &H) = 0;

        /// Objective gradient into the zeroed flat tangent vector @c grad.
        virtual Index eval_gradient(const Scalar *objective_scale, const Scalar *x,
                                    Scalar *grad) = 0;

        /// Inequality Jacobians, see @ref GraphJacobianView.
        virtual Index eval_jacobian_ineq(const Scalar *, GraphJacobianView &) { return 0; }

        /// Flat inequality-constraint values @c g_ineq(x).
        virtual Index eval_gineq(const Scalar *, Scalar *) { return 0; }

        /// Evaluate the objective.
        virtual Index eval_f(const Scalar *objective_scale, const Scalar *x, Scalar *res) = 0;

        /**
         * @brief Lower / upper bounds for the inequality constraints.
         *
         * Layout: per-block concatenation. Use @c +inf / @c -inf for
         * unbounded sides. Default: no-op (the @c GraphAbstract default
         * reports zero inequality constraints, so the bounds vectors are
         * empty).
         */
        virtual Index get_bounds(Scalar *, Scalar *) const { return 0; }

        /// Initial primal guess (flat vector, all blocks concatenated).
        virtual Index get_initial(Scalar *x) const = 0;

        /**
         * @brief Retraction @c x_next = retract(x, alpha * delta_x).
         *
         * For Euclidean blocks this is @c x_next = x + alpha * delta_x.
         *
         * @param x        Flat primal vector.
         * @param delta_x  Flat tangent vector.
         * @param alpha    Step size.
         * @param x_next   Output flat primal vector.
         */
        virtual void apply_retraction(const Scalar *x, const Scalar *delta_x, const Scalar alpha,
                                      Scalar *x_next) = 0;

        virtual ~GraphAbstract() = default;
    };
} // namespace fatrop

#endif // __fatrop_graph_graph_abstract_hpp__
