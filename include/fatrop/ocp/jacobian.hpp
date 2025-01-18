//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_jacobian_hpp__
#define __fatrop_ocp_jacobian_hpp__

#include "fatrop/nlp/jacobian.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/ocp/fwd.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include <vector>

/**
 * @file jacobian.hpp
 * @brief Defines the specialized Jacobian structure for Optimal Control Problems (OCPs).
 *
 * This file contains the specialization of the Jacobian structure template
 * for Optimal Control Problems (OCPs). It represents the constraint Jacobian
 * of the Karush-Kuhn-Tucker (KKT) system.
 *
 */

namespace fatrop
{
    typedef ProblemInfo<OcpType> OcpInfo;
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
        Jacobian(const ProblemDims<OcpType> &dims);

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
         *  BAbt[:nu, :] is reserved for control Jacobian blocks, while BAbt[nu:nu+nx, :] is
         * reserved for state Jacobian blocks. BAbt[-1, :] is reserved for the right-hand side.
         */
        std::vector<MatRealAllocated> BAbt;

        /**
         * @brief Constraint Jacobian of path equality constraints. Similar to BAbt, it also has
         * an additional row for the right-hand side.
         *
         * Matrix dimensions: (nu[k] + nx[k]) x ng[k]
         * Where:
         *   nx[k]: number of states at time step k
         *   nu[k]: number of control inputs at time step k
         *   ng[k]: number of equality constraints at time step k
         *  Gg_eqt[:nu, :] is reserved for control Jacobian blocks, while Gg_eqt[nu:nu+nx, :] is
         * reserved for state Jacobian blocks. Gg_eqt[-1, :] is reserved for the right-hand side.
         */
        std::vector<MatRealAllocated> Gg_eqt;

        /**
         * @brief Constraint Jacobian of path inequality constraints. Similar to BAbt, it also
         * has an additional row for the right-hand side.
         *
         * Matrix dimensions: (nu[k] + nx[k] + 1) x ng_ineq[k]
         * Where:
         *   nx[k]: number of states at time step k
         *   nu[k]: number of control inputs at time step k
         *   ng_ineq[k]: number of inequality constraints at time step k
         * Gg_ineqt[:nu, :] is reserved for control Jacobian blocks, while Gg_ineqt[nu:nu+nx, :] is
         * reserved for state Jacobian blocks. Gg_ineqt[-1, :] is reserved for the right-hand side.
         */
        std::vector<MatRealAllocated> Gg_ineqt;

        void apply_on_right(const OcpInfo& info, const VecRealView &x, Scalar alpha, const VecRealView& y, VecRealView &out) const;
        void transpose_apply_on_right(const OcpInfo& info, const VecRealView &mult_eq, Scalar alpha, const VecRealView& y, VecRealView &out) const;
        void get_rhs(const OcpInfo& info, VecRealView &rhs) const;
        void set_rhs(const OcpInfo& info, const VecRealView &rhs);
        // make printable 
        friend std::ostream &operator<<(std::ostream &os, const Jacobian &jac);
    };
} // namespace fatrop

#endif //__fatrop_ocp_jacobian_hpp__
