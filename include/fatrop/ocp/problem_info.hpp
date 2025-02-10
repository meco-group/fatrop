//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_problem_info__
#define __fatrop_ocp_problem_info__

#include "fatrop/nlp/problem_info.hpp"
#include "fatrop/ocp/dims.hpp"
#include "fatrop/ocp/type.hpp"
#include <vector>

namespace fatrop
{
    /**
     * @brief Specialization of ProblemInfo for Optimal Control Problems (OCP).
     *
     * This struct contains information about how variables and constraints of an OCP are ordered
     * internally.    */
    template <> struct ProblemInfo<OcpType>
    {
        /**
         * @brief Construct a new ProblemInfo object
         *
         * @param dims The dimensions of the OCP
         */
        ProblemInfo(const ProblemDims<OcpType> &dims);

        /// The dimensions of the OCP
        const ProblemDims<OcpType> dims;

        // Number of primal variables
        Index number_of_primal_variables;

        // Number of slack variables
        Index number_of_slack_variables;

        // Number of equality constraints
        Index number_of_eq_constraints; 

        /// Offset for the primal variables in the overall variable vector
        Index offset_primal;

        /**
         * @brief Offsets for the primal variables (u and x) at each stage
         *
         * The primal variables (x) are concatenated as x = [u_0, x_0, ..., u_{K-1}, x_{K-1}]
         */
        std::vector<Index> offsets_primal_u;
        std::vector<Index> offsets_primal_x;

        /// Offset for the slack variables in the overall variable vector (combined x and s vector)
        Index offset_slack;

        /**
         * @brief Offsets for the slack variables at each stage
         *
         * The slack variables (s) are concatenated as s = [s_0, s_1, ..., s_{K-1}],
         * where s_i is the slack variable related to the inequality path constraint at stage i.
         * With these offsets, the slack variables, the inequality constraints and the dual bound
         * multipliers are related.
         */
        std::vector<Index> offsets_slack; 
        std::vector<Index> offsets_eq;
        std::vector<Index> offsets_dyn;

        /// Offsets for the quantities related to equality constraints, the constraints are ordered
        /// as a concatenation of [g_eq_dyn, g_eq_path, g_eq_slack]

        /// Number of dynamics equality constraints
        Index number_of_g_eq_dyn;
        /// Number of path equality constraints
        Index number_of_g_eq_path;
        /// Number of slack equality constraints
        Index number_of_g_eq_slack;
        /// Offset for the dynamics equality constraints
        Index offset_g_eq_dyn;

        /// Offset for the path equality constraints
        Index offset_g_eq_path;

        /// Offset for the slack equality constraints
        Index offset_g_eq_slack;

        /**
         * @brief Offsets for the dynamics equality constraints at each stage
         *
         * The equality dynamics constraints are concatenated as
         * g_eq_dyn = [x_1 - f(x_0, u_0), x_2 - f(x_1, u_1), ..., x_{K-1} - f(x_{K-2}, u_{K-2})]
         * where f is the dynamics function. With these offsets, the dynamics equality constraints
         * and the dual equality multipliers are related.
         */
        std::vector<Index> offsets_g_eq_dyn;

        /**
         * @brief Offsets for the path equality constraints at each stage
         *
         * The equality path constraints are concatenated as
         * g_eq_path = [g_0(x_0, u_0), g_1(x_1, u_1), ..., g_{K-1}(x_{K-1}, u_{K-1})].
         * With these offsets, the path equality constraints and the dual equality multipliers are
         * related.
         */
        std::vector<Index> offsets_g_eq_path;

        /**
         * @brief Offsets for the slack equality constraints at each stage
         *
         * The slack equality constraints are concatenated as
         * g_eq_slack = [g_ineq_0(x_0, u_0) - s_0, g_ineq_1(x_1, u_1) - s_1, ...,
         * g_ineq_{K-1}(x_{K-1}, u_{K-1}) - s_{K-1}].
         * With these offsets, the slack equality constraints and the dual equality multipliers are
         * related.
         *
         * These constraints are added to rewrite L<=g_ineq(x)<=U equivalently as g_ineq(x) - s = 0
         * and L<=s<=U.
         */
        std::vector<Index> offsets_g_eq_slack;

        std::vector<Index> number_of_stage_variables;

        Index pd_orig_offset_primal;
        Index pd_orig_offset_slack;
        Index pd_orig_offset_mult;
        Index pd_orig_offset_zl;
        Index pd_orig_offset_zu;

        // Info for the restoration phase
        Index number_of_slack_variables_resto;

        Index pd_resto_offset_primal;
        Index pd_resto_offset_slack;
        Index pd_resto_offset_mult;
        Index pd_resto_offset_zl;
        Index pd_resto_offset_zu;

        Index pd_resto_offset_zp;
        Index pd_resto_offset_zn;
        // when normal and slack variables are in the same vector
        Index offset_slack_p;
        Index offset_slack_n;
        // Offset when in vector with all slack variables
        Index offset_s;
        Index offset_p;
        Index offset_n;
    };
} // namespace fatrop

#endif // __fatrop_ocp_problem_info__
