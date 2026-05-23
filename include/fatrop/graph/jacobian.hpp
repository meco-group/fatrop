//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_jacobian_hpp__
#define __fatrop_graph_jacobian_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/graph/dims.hpp"
#include "fatrop/graph/fwd.hpp"
#include "fatrop/graph/problem_type.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/nlp/jacobian.hpp"

#include <iosfwd>
#include <vector>

namespace fatrop
{
    typedef ProblemInfo<GraphProblem> GraphInfo;

    /**
     * @brief Inequality-constraint Jacobian specialisation for graph problems.
     *
     * The inequality constraints are block-local: each row depends only on the
     * variables of exactly one block @c k. The Jacobian is therefore stored
     * per block as a transposed @c (n_k + 1) × ng_ineq[k] BLASFEO matrix —
     * the top @c n_k rows hold the Jacobian columns (in tangent space), the
     * trailing row holds the right-hand side @c g_ineq_k(x_k). This is the
     * same layout convention used by the OCP and dense Jacobian
     * specialisations.
     *
     * Graph problems have no equality constraints, so no @c Gg_eqt vector is
     * needed.
     */
    template <> struct Jacobian<GraphProblem>
    {
        Jacobian(const ProblemDims<GraphProblem> &dims);

        /// Per-block inequality Jacobian (transposed + RHS row).
        /// @c Gg_ineqt[k] has shape @c (n_k + 1) × ng_ineq[k]. If
        /// @c ng_ineq[k] == 0 the entry is a zero-column matrix.
        std::vector<MatRealAllocated> Gg_ineqt;

        void apply_on_right(const GraphInfo &info, const VecRealView &x, Scalar alpha,
                            const VecRealView &y, VecRealView &out) const;
        void transpose_apply_on_right(const GraphInfo &info, const VecRealView &mult_eq,
                                      Scalar alpha, const VecRealView &y, VecRealView &out) const;
        void get_rhs(const GraphInfo &info, VecRealView &rhs) const;
        void set_rhs(const GraphInfo &info, const VecRealView &rhs);

        friend std::ostream &operator<<(std::ostream &os, const Jacobian<GraphProblem> &jac);
    };
} // namespace fatrop

#endif // __fatrop_graph_jacobian_hpp__
