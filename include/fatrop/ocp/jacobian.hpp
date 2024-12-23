//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_jacobian_hpp__
#define __fatrop_ocp_jacobian_hpp__

#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/nlp/jacobian.hpp"
#include "fatrop/ocp/dims.hpp"
#include "fatrop/ocp/type.hpp"

/**
 * @file jacobian.hpp
 * @brief Defines the specialized Jacobian structure for Optimal Control Problems (OCPs).
 *
 * This file contains the specialization of the Jacobian structure template
 * for Optimal Control Problems (OCPs). It represents the constraint Jacobian
 * of the Karush-Kuhn-Tucker (KKT) system, which is crucial for solving OCPs.
 *
 * Optimal Control Problems involve finding a control strategy that optimizes
 * a performance criterion subject to dynamic constraints and other path constraints.
 * The Jacobian structure defined here is tailored to efficiently handle the
 * specific sparsity patterns and block structures that arise in OCPs.
 */

namespace fatrop
{
    /**
     * @brief Specialization of the Jacobian structure for Optimal Control Problems.
     *
     * This specialization of the Jacobian structure is designed specifically
     * for Optimal Control Problems (OCPs). It represents the constraint Jacobian
     * of the Karush-Kuhn-Tucker (KKT) system for OCPs, taking into account their
     * specific structure and requirements.
     *
     * The Jacobian is structured to efficiently handle the block-sparse nature
     * of OCPs, which typically involve dynamics constraints, path constraints,
     * and terminal constraints over multiple time steps.
     */
    template <> struct Jacobian<OcpType>
    {
        /**
         * @brief Construct a new Jacobian object for an OCP.
         *
         * @param dims The dimensions of the OCP, used to allocate appropriate memory.
         */
        Jacobian(const OcpDims &dims);

        /**
         * @brief Constraint Jacobian of the dynamics.
         *
         * The Jacobian is represented by a transposed matrix with an additional row
         * for the right-hand side. This structure allows for efficient simultaneous
         * factorization and solve in the Riccati recursion.
         *
         * Matrix dimensions: (nu[k] + nx[k] + 1) x nx[k+1]
         * Where:
         *   nu[k]: number of control inputs at time step k
         *   nx[k]: number of states at time step k
         *   nx[k+1]: number of states at time step k+1
         *  jac_dyn[:nu, :] is reserved for control Jacobian blocks, while jac_dyn[nu:nu+nx, :] is
         * reserved for state Jacobian blocks. jac_dyn[-1, :] is reserved for the right-hand side.
         */
        std::vector<MatRealAllocated> jac_dyn;

        /**
         * @brief Constraint Jacobian of path equality constraints. Similar to jac_dyn, it also has
         * an additional row for the right-hand side.
         *
         * Matrix dimensions: (nu[k] + nx[k]) x ng[k]
         * Where:
         *   nx[k]: number of states at time step k
         *   nu[k]: number of control inputs at time step k
         *   ng[k]: number of equality constraints at time step k
         *  jac_eq[:nu, :] is reserved for control Jacobian blocks, while jac_eq[nu:nu+nx, :] is
         * reserved for state Jacobian blocks. jac_eq[-1, :] is reserved for the right-hand side.
         */
        std::vector<MatRealAllocated> jac_eq;

        /**
         * @brief Constraint Jacobian of path inequality constraints. Similar to jac_dyn, it also
         * has an additional row for the right-hand side.
         *
         * Matrix dimensions: (nu[k] + nx[k] + 1) x ng_ineq[k]
         * Where:
         *   nx[k]: number of states at time step k
         *   nu[k]: number of control inputs at time step k
         *   ng_ineq[k]: number of inequality constraints at time step k
         * jac_ineq[:nu, :] is reserved for control Jacobian blocks, while jac_ineq[nu:nu+nx, :] is
         * reserved for state Jacobian blocks. jac_ineq[-1, :] is reserved for the right-hand side.
         */
        std::vector<MatRealAllocated> jac_ineq;
    };
} // namespace fatrop

#endif //__fatrop_ocp_jacobian_hpp__
