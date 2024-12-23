//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_hessian_hpp__
#define __fatrop_ocp_hessian_hpp__

#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/nlp/hessian.hpp"
#include "fatrop/ocp/dims.hpp"
#include "fatrop/ocp/type.hpp"

/**
 * @file hessian.hpp
 * @brief Defines the specialized Hessian structure for Optimal Control Problems (OCPs).
 *
 * This file contains the specialization of the hessian structure template
 * for Optimal Control Problems (OCPs). It represents the constraint Hessian
 * of the Karush-Kuhn-Tucker (KKT) system, which is crucial for solving OCPs.
 *
 * Optimal Control Problems involve finding a control strategy that optimizes
 * a performance criterion subject to dynamic constraints and other path constraints.
 * The Hessian structure defined here is tailored to efficiently handle the
 * specific sparsity patterns and block structures that arise in OCPs.
 */

namespace fatrop
{
    /**
     * @brief Specialization of the Hessian structure for Optimal Control Problems.
     *
     * This specialization of the Hessian structure is designed specifically
     * for Optimal Control Problems (OCPs). It represents the Hessian
     * of the Karush-Kuhn-Tucker (KKT) system for OCPs, taking into account their
     * specific structure and requirements.
     *
     * The Hessian is structured to efficiently handle the block-sparse nature
     * of OCPs, which typically involve dynamics constraints, path constraints,
     * and terminal constraints over multiple time steps.
     */
    template <> struct Hessian<OcpType>
    {
        /**
         * @brief Construct a new Hessian object for an OCP.
         *
         * @param dims The dimensions of the OCP, used to allocate appropriate memory.
         */
        Hessian(const OcpDims &dims);

        /**
         * @brief Constraint Hessian of the dynamics.
         *
         * The Hessian is represented by a transposed matrix with an additional row
         * for the right-hand side. This structure allows for efficient simultaneous
         * factorization and solve in the Riccati recursion.
         *
         * Matrix dimensions: (nu[k] + nx[k] + 1) x (nu[k] + nx[k] + 1)
         * Where:
         *   nu[k]: number of control inputs at time step k
         *   nx[k]: number of states at time step k
         *   nx[k+1]: number of states at time step k+1
         *   hess[:nu, :nu] is reserved for control-control Hessian blocks,
         *   hess[nu:nu+nx, nu:nu+nx] is reserved for state-state Hessian blocks,
         *   hess[:nu, nu:nu+nx] = hess[nu:nu+nx, :nu]^T is reserved for control-state "skew"
         *   hess[-1, :] is reserved for the right-hand side.
         */
        std::vector<MatRealAllocated> hess;
    };
} // namespace fatrop

#endif //__fatrop_ocp_jacobian_hpp__
