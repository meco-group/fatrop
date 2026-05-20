//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_dims_hpp__
#define __fatrop_ocp_dims_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/nlp/dims.hpp"
#include "fatrop/ocp/type.hpp"
#include <vector>

namespace fatrop
{
    /**
     * @brief Represents the dimensions of an Optimal Control Problem (OCP).
     *
     * This structure holds the dimensions of various components in an OCP,
     * including the number of states, controls, and constraints for each stage.
     */
    template <>
    struct ProblemDims<OcpType>
    {
        /**
         * @brief Constructs an ProblemDims object (lvalue reference version).
         *
         * @param K Number of stages / time steps in the Ocp. This means that there are K state
         * vectors, and K-1 control vectors used in the dynamics.
         * @param nu Vector of number of controls for each stage.
         * @param nx Vector of number of states for each stage.
         * @param ng Vector of number of equality constraints for each stage.
         * @param ng_ineq Vector of number of inequality constraints for each stage.
         */
        ProblemDims(int K, const std::vector<Index> &nu, const std::vector<Index> &nx,
                const std::vector<Index> &ng, const std::vector<Index> &ng_ineq);

        /**
         * @brief Constructs an ProblemDims object (rvalue reference version).
         *
         * @param k Number of stages in the OCP.
         * @param nu Vector of number of controls for each stage.
         * @param nx Vector of number of states for each stage.
         * @param ng Vector of number of equality constraints for each stage.
         * @param ng_ineq Vector of number of inequality constraints for each stage.
         */
        ProblemDims(int K, std::vector<Index> &&nu, std::vector<Index> &&nx, std::vector<Index> &&ng,
                std::vector<Index> &&ng_ineq);
        /**
         * @brief Constructs a ProblemDims object with separate tangent-space dimensions.
         *
         * Use this overload for Lie-group / manifold optimization where the primal state /
         * control lives on a manifold of dimension @c nx[k] / @c nu[k] but the
         * search-direction (tangent) lives in a space of dimension @c nx_tan[k] /
         * @c nu_tan[k]. Pass equal vectors (or use the other constructors) for standard
         * Euclidean problems.
         */
        ProblemDims(int K, const std::vector<Index> &nu, const std::vector<Index> &nx,
                    const std::vector<Index> &ng, const std::vector<Index> &ng_ineq,
                    const std::vector<Index> &nu_tan, const std::vector<Index> &nx_tan);
        void check_problem_dimensions() const;

        const int K;                                  ///< Number of stages in the OCP.
        const std::vector<Index> number_of_controls; ///< Number of controls for each stage.
        const std::vector<Index> number_of_states;    ///< Number of states for each stage.
        const std::vector<Index>
            number_of_eq_constraints; ///< Number of path equality constraints for each stage.
        const std::vector<Index>
            number_of_ineq_constraints; ///< Number of path inequality constraints for each stage.
        /// Tangent-space dimension of the controls at each stage. Equal to
        /// @c number_of_controls unless an OCP overrides @c get_nu_tangent.
        const std::vector<Index> number_of_tangent_controls;
        /// Tangent-space dimension of the state at each stage. Equal to
        /// @c number_of_states unless an OCP overrides @c get_nx_tangent.
        const std::vector<Index> number_of_tangent_states;
    };

} // namespace fatrop

#endif //__fatrop_ocp_dims_hpp__
