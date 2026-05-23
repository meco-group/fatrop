//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_problem_type_hpp__
#define __fatrop_graph_problem_type_hpp__

/**
 * @file problem_type.hpp
 * @brief Defines the GraphProblem tag for block-structured NLPs.
 *
 * A graph problem has the form
 *
 *     min   f(x)
 *     s.t.  L_k <= g_k_i(x_k) <= U_k_i   for each block k, constraint i
 *
 * where the primal variable @c x is partitioned into @c N blocks
 * @c x_0, ..., x_{N-1}, the Lagrangian Hessian has a user-supplied symmetric
 * block-sparse structure (described by an undirected adjacency graph), and the
 * inequality constraints are block-local: each constraint depends on the
 * variables of exactly one block. The problem has no equality constraints.
 *
 * The corresponding KKT system, after the standard Schur elimination of the
 * slack variables and bound multipliers, reduces to an SPD block-sparse system
 * with the same block sparsity as the Hessian. It is solved by
 * @ref BlockCholeskySolver.
 *
 * @note @c GraphProblem is the problem-class tag for the interior-point
 *       algorithm (analogue of @c OcpType and @c DenseType); it is distinct
 *       from @ref GraphType, which tags only the SPD block linear solver.
 */

namespace fatrop
{
    struct GraphProblem
    {
    };
} // namespace fatrop

#endif // __fatrop_graph_problem_type_hpp__
