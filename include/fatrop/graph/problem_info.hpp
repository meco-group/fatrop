//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_problem_info_hpp__
#define __fatrop_graph_problem_info_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/graph/block_sparsity.hpp"
#include "fatrop/graph/dims.hpp"
#include "fatrop/graph/problem_type.hpp"
#include "fatrop/nlp/problem_info.hpp"

#include <memory>
#include <vector>

namespace fatrop
{
    /**
     * @brief Specialisation of @c ProblemInfo for graph-structured NLPs.
     *
     * The layout mirrors the dense / OCP variants: a single primal block (in
     * tangent space, equal to the primal space for Euclidean problems) followed
     * by a slack block, an equality-multiplier block (empty for graph
     * problems), and lower / upper bound-multiplier blocks. Only the fields
     * actually consumed by the (problem-type-agnostic) IP-algorithm code and by
     * the graph PD-system code are stored.
     */
    template <> struct ProblemInfo<GraphProblem>
    {
        ProblemInfo(const ProblemDims<GraphProblem> &dims);

        const ProblemDims<GraphProblem> dims;

        /// Block sparsity pattern of the Hessian. Owned because it is needed
        /// by the linear solver after construction; references would dangle if
        /// the original @c dims object went out of scope first.
        std::shared_ptr<BlockSparsityPattern> sparsity;

        Index number_of_primal_variables;   ///< == dims.nx
        Index number_of_tangent_variables;  ///< == dims.nx (no manifold variables yet)

        Index number_of_slack_variables;    ///< == dims.ng_ineq_total
        Index offset_slack;                 ///< offset of the slack damping inside D_x

        Index number_of_eq_constraints;     ///< == dims.ng_ineq_total (only slack-eq constraints)
        Index number_of_g_eq_slack;         ///< == dims.ng_ineq_total
        Index offset_g_eq_slack;            ///< 0 — no path-eq block ahead

        // PD-system (original) block offsets: layout is
        //   [tangent | slack | mult | zl | zu].
        Index pd_orig_offset_primal;
        Index pd_orig_offset_slack;
        Index pd_orig_offset_mult;
        Index pd_orig_offset_zl;
        Index pd_orig_offset_zu;

        // Restoration phase: graph problems have no restoration phase, but the
        // shared IP-algorithm code reads these fields, so we initialise them to
        // safe defaults that never get exercised.
        Index number_of_slack_variables_resto;
        Index pd_resto_offset_primal;
        Index pd_resto_offset_slack;
        Index pd_resto_offset_mult;
        Index pd_resto_offset_zl;
        Index pd_resto_offset_zu;
        Index pd_resto_offset_zp;
        Index pd_resto_offset_zn;

        Index offset_slack_p;
        Index offset_slack_n;
        Index offset_s;
        Index offset_p;
        Index offset_n;

        std::vector<bool> constraint_allows_dual_damping;
    };
} // namespace fatrop

#endif // __fatrop_graph_problem_info_hpp__
