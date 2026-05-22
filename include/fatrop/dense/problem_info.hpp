//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_dense_problem_info_hpp__
#define __fatrop_dense_problem_info_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/dense/dims.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/nlp/problem_info.hpp"
#include <vector>

namespace fatrop
{
    /**
     * @brief Specialization of ProblemInfo for dense NLP problems.
     *
     * A dense problem has a single primal vector and no time / staging structure,
     * so most of the OCP-flavored offsets (controls, tangent space, dynamics block)
     * collapse to trivial values. We keep only what is actually consumed by the
     * (problem-type agnostic) IP-algorithm code and by the dense PD system code:
     *
     *  - the variable / slack / equality counts,
     *  - the equality-slack block offset (used by the @c constr_viol_ineq
     *    specialization),
     *  - the slack offset inside the combined (primal, slack) damping vector,
     *  - the original- and restoration-phase PD-system block offsets,
     *  - the restoration-phase @c p / @c n slack bookkeeping.
     *
     * Anything implied by the dense layout (primal/tangent/path-eq all live at
     * offset 0, no dynamics block) is not stored — the dense code uses 0 directly.
     */
    template <> struct ProblemInfo<DenseType>
    {
        ProblemInfo(const ProblemDims<DenseType> &dims);

        const ProblemDims<DenseType> dims;

        Index number_of_primal_variables;  ///< == dims.nx
        Index number_of_tangent_variables; ///< == dims.nx (no manifold variables)

        Index number_of_slack_variables; ///< == dims.ng_ineq
        Index offset_slack;              ///< offset of the slack damping inside D_x

        Index number_of_eq_constraints; ///< ng + ng_ineq
        Index number_of_g_eq_slack;     ///< == dims.ng_ineq
        Index offset_g_eq_slack;        ///< == dims.ng (offset of slack-eq block in eq vector)

        // PD-system (original) block offsets.
        Index pd_orig_offset_primal;
        Index pd_orig_offset_slack;
        Index pd_orig_offset_mult;
        Index pd_orig_offset_zl;
        Index pd_orig_offset_zu;

        // Restoration phase.
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

#endif // __fatrop_dense_problem_info_hpp__
