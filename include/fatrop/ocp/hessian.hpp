//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_hessian_hpp__
#define __fatrop_ocp_hessian_hpp__

#include "fatrop/nlp/hessian.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/ocp/fwd.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include <vector>

/**
 * @file hessian.hpp
 * @brief Defines the specialized Hessian structure for Optimal Control Problems (OCPs).
 *
 * This file contains the specialization of the hessian structure template
 * for Optimal Control Problems (OCPs). It represents the constraint Hessian
 * of the Karush-Kuhn-Tucker (KKT) system.
 *
 */

namespace fatrop
{
    typedef ProblemInfo<OcpType> OcpInfo;
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
         *   RSQrqt[:nu, :nu] is reserved for control-control Hessian blocks,
         *   RSQrqt[nu:nu+nx, nu:nu+nx] is reserved for state-state Hessian blocks,
         *   RSQrqt[:nu, nu:nu+nx] = RSQrqt[nu:nu+nx, :nu]^T is reserved for control-state "skew"
         *   RSQrqt[-1, :] is reserved for the right-hand side.
         */
        std::vector<MatRealAllocated> RSQrqt;

        void apply_on_right(const OcpInfo& info, const VecRealView& x, Scalar alpha, const VecRealView& y, VecRealView& out) const;
        void get_rhs(const OcpInfo& info, VecRealView& out) const;
        void set_rhs(const OcpInfo& info, const VecRealView& in);
    };
} // namespace fatrop

#endif //__fatrop_ocp_jacobian_hpp__
