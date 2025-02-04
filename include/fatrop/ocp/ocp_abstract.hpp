//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_abstract_hpp__
#define __fatrop_ocp_abstract_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/ocp/ocp_eval.hpp"

/**
 * @file ocp_abstract.hpp
 * @brief Defines the abstract base class for Optimal Control Problems (OCPs).
 *
 * This file contains the abstract base class OcpAbstract, which provides
 * the interface for defining and solving Optimal Control Problems (OCPs).
 * It includes methods for evaluating dynamics, constraints, objective functions,
 * and their derivatives, as well as methods for retrieving problem dimensions
 * and initial guesses.
 */

namespace fatrop
{
    // tag for OcpAbstract
    struct OcpAbstractDynamic
    {
    };
    /**
     * @brief Abstract base class for Optimal Control Problems.
     *
     * This class defines the interface for Optimal Control Problems (OCPs).
     * It provides methods for evaluating dynamics, constraints, objective functions,
     * and their derivatives. It also includes methods for retrieving problem
     * dimensions, initial guesses, and bounds.
     *
     * Derived classes should implement these methods to define specific OCPs.
     */

    template <> class OcpAbstractTpl<OcpAbstractDynamic>
    {
    public:
        /**
         * @brief Get the number of states for a given time step.
         * @param k Time step
         * @return Number of states
         */
        virtual Index get_nx(const Index k) const = 0;

        /**
         * @brief Get the number of inputs for a given time step.
         * @param k Time step
         * @return Number of inputs
         */
        virtual Index get_nu(const Index k) const = 0;

        /**
         * @brief Get the number of equality constraints for a given time step.
         * @param k Time step
         * @return Number of equality constraints
         */
        virtual Index get_ng(const Index k) const = 0;

        /**
         * @brief Get the number of inequality constraints for a given time step.
         * @param k Time step
         * @return Number of inequality constraints
         */
        virtual Index get_ng_ineq(const Index k) const = 0;

        /**
         * @brief Get the horizon length of the OCP.
         * @return Horizon length
         */
        virtual Index get_horizon_length() const = 0;
        /**
         * @brief Evaluate the discretized dynamics.
         *
         * This method evaluates the vertical concatenation of A_k^T and B_k^T
         * from the linearized dynamics x_{k+1} = A_k x_k + B_k u_k + b_k.
         * The matrix is a BLASFEO matrix type.
         *
         * Note: The last row of the matrix (b_k^T) is used internally and should not
         * be evaluated by the user.
         *
         * @param states_kp1 Pointer to nx_{k+1}-array states of time step k+1
         * @param inputs_k Pointer to array inputs of time step k
         * @param states_k Pointer to array states of time step k
         * @param res Pointer to (nu+nx+1 x nu+nx)-matrix
         * @param k Time step
         * @return Status code
         */
        virtual Index eval_BAbt(const Scalar *states_kp1, const Scalar *inputs_k,
                                const Scalar *states_k, MAT *res, const Index k) = 0;

        /**
         * @brief Evaluate the stagewise Lagrangian Hessian.
         *
         * This method evaluates the vertical concatenation of:
         * 1. The Hessian of the Lagrangian to the concatenation of (u_k, x_k)
         * 2. The first order derivative of the Lagrangian Hessian to the concatenation of (u_k,
         * x_k) The matrix is a BLASFEO matrix type.
         *
         * Note: The last row of the matrix is used internally and should not
         * be evaluated by the user.
         *
         * @param objective_scale Scale factor for objective function (usually 1.0)
         * @param inputs_k Pointer to array inputs of time step k
         * @param states_k Pointer to array states of time step k
         * @param lam_dyn_k Pointer to array dual variables for dynamics of time step k
         * @param lam_eq_k Pointer to array dual variables for equality constraints of time step k
         * @param lam_eq_ineq_k Pointer to array dual variables for inequality constraints of time
         * step k
         * @param res Pointer to (nu+nx+1 x nu+nx)-matrix
         * @param k Time step
         * @return Status code
         */
        virtual Index eval_RSQrqt(const Scalar *objective_scale, const Scalar *inputs_k,
                                  const Scalar *states_k, const Scalar *lam_dyn_k,
                                  const Scalar *lam_eq_k, const Scalar *lam_eq_ineq_k, MAT *res,
                                  const Index k) = 0;

        /**
         * @brief Evaluate the stagewise equality constraints Jacobian.
         *
         * This method evaluates the Jacobian of the equality constraints
         * with respect to the concatenation of (u_k, x_k).
         * The matrix is a BLASFEO matrix type.
         *
         * Note: The last row of the matrix is used internally and should not
         * be evaluated by the user.
         *
         * @param inputs_k Pointer to array inputs of time step k
         * @param states_k Pointer to array states of time step k
         * @param res Pointer to (nu+nx+1 x ng)-matrix
         * @param k Time step
         * @return Status code
         */
        virtual Index eval_Ggt(const Scalar *inputs_k, const Scalar *states_k, MAT *res,
                               const Index k) = 0;

        /**
         * @brief Evaluate the stagewise inequality constraints Jacobian.
         *
         * This method evaluates the Jacobian of the inequality constraints
         * with respect to the concatenation of (u_k, x_k).
         * The matrix is a BLASFEO matrix type.
         *
         * Note: The last row of the matrix is used internally and should not
         * be evaluated by the user.
         *
         * @param inputs_k Pointer to array inputs of time step k
         * @param states_k Pointer to array states of time step k
         * @param res Pointer to (nu+nx+1 x ng_ineq)-matrix, BLASFEO matrix type
         * @param k Time step
         * @return Status code
         */
        virtual Index eval_Ggt_ineq(const Scalar *inputs_k, const Scalar *states_k, MAT *res,
                                    const Index k) = 0;
        /**
         * @brief Evaluate the dynamics constraint violation.
         *
         * This method evaluates b_k = -x_{k+1} + f_k(u_k, x_k, p_k, p).
         *
         * @param states_kp1 Pointer to array states of time step k+1
         * @param inputs_k Pointer to array inputs of time step k
         * @param states_k Pointer to array states of time step k
         * @param res Pointer to array nx_{k+1}-vector
         * @param k Time step
         * @return Status code
         */
        virtual Index eval_b(const Scalar *states_kp1, const Scalar *inputs_k,
                             const Scalar *states_k, Scalar *res, const Index k) = 0;

        /**
         * @brief Evaluate the equality constraint violation.
         *
         * This method evaluates g_k = g_k(u_k, x_k, p_k, p).
         *
         * @param inputs_k Pointer to array inputs of time step k
         * @param states_k Pointer to array states of time step k
         * @param res Pointer to array ng-vector
         * @param k Time step
         * @return Status code
         */
        virtual Index eval_g(const Scalar *inputs_k, const Scalar *states_k, Scalar *res,
                             const Index k) = 0;

        /**
         * @brief Evaluate the inequality constraint violation.
         *
         * This method evaluates g_ineq_k = g_ineq_k(u_k, x_k, p_k, p).
         *
         * @param inputs_k Pointer to array inputs of time step k
         * @param states_k Pointer to array states of time step k
         * @param res Pointer to array ng_ineq-vector
         * @param k Time step
         * @return Status code
         */
        virtual Index eval_gineq(const Scalar *inputs_k, const Scalar *states_k, Scalar *res,
                                 const Index k) = 0;

        /**
         * @brief Evaluate the gradient of the objective function.
         *
         * This method evaluates the gradient of the objective function (not the Lagrangian!)
         * with respect to the concatenation of (u_k, x_k).
         *
         * @param objective_scale Pointer to objective scale
         * @param inputs_k Pointer to array inputs of time step k
         * @param states_k Pointer to array states of time step k
         * @param res Pointer to (nu+nx)-array
         * @param k Time step
         * @return Status code
         */
        virtual Index eval_rq(const Scalar *objective_scale, const Scalar *inputs_k,
                              const Scalar *states_k, Scalar *res, const Index k) = 0;

        /**
         * @brief Evaluate the objective function value.
         *
         * @param objective_scale Pointer to array objective scale
         * @param inputs_k Pointer to array inputs of time step k
         * @param states_k Pointer to array states of time step k
         * @param res Pointer to Scalar
         * @param k Time step
         * @return Status code
         */
        virtual Index eval_L(const Scalar *objective_scale, const Scalar *inputs_k,
                             const Scalar *states_k, Scalar *res, const Index k) = 0;

        /**
         * @brief Get the bounds of the inequalities at a given stage.
         *
         * @param lower Pointer to ng_ineq-vector for lower bounds
         * @param upper Pointer to ng_ineq-vector for upper bounds
         * @param k Time step
         * @return Status code
         */
        virtual Index get_bounds(Scalar *lower, Scalar *upper, const Index k) const = 0;

        /**
         * @brief Get the default initial guess for the states of a given stage.
         *
         * @param xk Pointer to states of time step k
         * @param k Time step
         * @return Status code
         */
        virtual Index get_initial_xk(Scalar *xk, const Index k) const = 0;

        /**
         * @brief Get the default initial guess for the inputs of a given stage.
         *
         * @param uk Pointer to inputs of time step k
         * @param k Time step
         * @return Status code
         */
        virtual Index get_initial_uk(Scalar *uk, const Index k) const = 0;

        /**
         * @brief Virtual destructor for OcpAbstract.
         *
         * This virtual destructor ensures proper cleanup of derived classes.
         */
        virtual ~OcpAbstractTpl() = default;
    };

    typedef OcpAbstractTpl<OcpAbstractDynamic> OcpAbstract;

} // namespace fatrop

#endif //__fatrop_ocp_abstract_hpp__
