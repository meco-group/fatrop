//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_hessian_hpp__
#define __fatrop_graph_hessian_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/graph/block_pd_matrix.hpp"
#include "fatrop/graph/dims.hpp"
#include "fatrop/graph/fwd.hpp"
#include "fatrop/graph/problem_type.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/nlp/hessian.hpp"

#include <iosfwd>

namespace fatrop
{
    typedef ProblemInfo<GraphProblem> GraphInfo;

    /**
     * @brief Lagrangian Hessian specialisation for graph-structured problems.
     *
     * Stores:
     *   - @c H: the block-sparse symmetric Lagrangian Hessian (a @ref
     *     BlockPdMatrix, sharing the sparsity owned by @c ProblemDims);
     *   - @c grad: the dense objective gradient in tangent space, length
     *     @c dims.nx (concatenation of per-block gradients in the standard
     *     block-offset layout).
     *
     * Splitting @c H and @c grad mirrors the role of the @c (nx + 1) x nx
     * @c Hht matrix used by the dense / OCP variants: the top block is the
     * Hessian, the trailing row stores the gradient. Here the Hessian is
     * block-sparse so the gradient is kept as a separate vector instead of a
     * trailing row.
     */
    template <> struct Hessian<GraphProblem>
    {
        Hessian(const ProblemDims<GraphProblem> &dims);

        /// Block-sparse symmetric Lagrangian Hessian, in the user-supplied
        /// block sparsity pattern (lower triangle).
        BlockPdMatrix H;

        /// Objective gradient (length @c dims.nx).
        VecRealAllocated grad;

        void apply_on_right(const GraphInfo &info, const VecRealView &x, Scalar alpha,
                            const VecRealView &y, VecRealView &out) const;
        void get_rhs(const GraphInfo &info, VecRealView &out) const;
        void set_rhs(const GraphInfo &info, const VecRealView &in);
        void set_zero();

        friend std::ostream &operator<<(std::ostream &os, const Hessian<GraphProblem> &hess);
    };
} // namespace fatrop

#endif // __fatrop_graph_hessian_hpp__
