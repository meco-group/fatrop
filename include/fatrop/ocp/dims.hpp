//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_dims_hpp__
#define __fatrop_ocp_dims_hpp__

#include "fatrop/context/context.hpp"
#include <vector>

namespace fatrop
{
    /**
     * @brief Represents the dimensions of an Optimal Control Problem (OCP).
     *
     * This structure holds the dimensions of various components in an OCP,
     * including the number of states, controls, and constraints for each stage.
     */
    struct OcpDims
    {
        /**
         * @brief Constructs an OcpDims object (lvalue reference version).
         *
         * @param K Number of stages / time steps in the Ocp. This means that there are K state
         * vectors, and K-1 control vectors used in the dynamics.
         * @param nu Vector of number of controls for each stage.
         * @param nx Vector of number of states for each stage.
         * @param ng Vector of number of equality constraints for each stage.
         * @param ng_ineq Vector of number of inequality constraints for each stage.
         */
        OcpDims(int K, const std::vector<Index> &nu, const std::vector<Index> &nx,
                const std::vector<Index> &ng, const std::vector<Index> &ng_ineq);

        /**
         * @brief Constructs an OcpDims object (rvalue reference version).
         *
         * @param k Number of stages in the OCP.
         * @param nu Vector of number of controls for each stage.
         * @param nx Vector of number of states for each stage.
         * @param ng Vector of number of equality constraints for each stage.
         * @param ng_ineq Vector of number of inequality constraints for each stage.
         */
        OcpDims(int K, std::vector<Index> &&nu, std::vector<Index> &&nx, std::vector<Index> &&ng,
                std::vector<Index> &&ng_ineq);
        void check_problem_dimensions() const;

        const int K;                                  ///< Number of stages in the OCP.
        const std::vector<Index> number_of_controls; ///< Number of controls for each stage.
        const std::vector<Index> number_of_states;    ///< Number of states for each stage.
        const std::vector<Index>
            number_of_eq_constraints; ///< Number of path equality constraints for each stage.
        const std::vector<Index>
            number_of_ineq_constraints; ///< Number of path inequality constraints for each stage.
    };

} // namespace fatrop

#endif //__fatrop_ocp_dims_hpp__
