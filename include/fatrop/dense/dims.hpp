//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_dense_dims_hpp__
#define __fatrop_dense_dims_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/nlp/dims.hpp"

namespace fatrop
{
    /**
     * @brief Dimensions of a dense (unstructured) NLP problem.
     *
     * The primal variable @c x has size @c nx. The search direction lives in a
     * (possibly smaller) tangent space of size @c nx_tangent — defaults to
     * @c nx for plain Euclidean problems and is overridden for Lie-group
     * optimisation (e.g. unit quaternion with @c nx == 4, @c nx_tangent == 3).
     * There are @c ng equality and @c ng_ineq inequality constraints.
     */
    template <> struct ProblemDims<DenseType>
    {
        /// Euclidean problem: tangent dimension defaults to @c nx.
        ProblemDims(Index nx, Index ng, Index ng_ineq);

        /// Manifold / Lie-group problem: explicit tangent dimension.
        ProblemDims(Index nx, Index nx_tangent, Index ng, Index ng_ineq);

        void check_problem_dimensions() const;

        const Index nx;         ///< Size of the primal vector x.
        const Index nx_tangent; ///< Size of the search direction delta_x.
        const Index ng;         ///< Number of equality constraints.
        const Index ng_ineq;    ///< Number of inequality constraints.
    };

} // namespace fatrop

#endif // __fatrop_dense_dims_hpp__
