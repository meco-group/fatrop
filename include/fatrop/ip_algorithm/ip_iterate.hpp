//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_iterate_hpp__
#define __fatrop_ip_algorithm_ip_iterate_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/nlp/fwd.hpp"
#include "fatrop/nlp/hessian.hpp"
#include "fatrop/nlp/jacobian.hpp"
#include <memory>
#include <utility>
#include <vector>

namespace fatrop
{
    /**
     * @brief Represents an iterate in the Interior Point (IP) algorithm for a given problem type.
     *
     * This class encapsulates the state and operations related to a single iteration
     * in the IP algorithm, including primal and dual variables, search directions,
     * and various computed quantities.
     *
     * @tparam ProblemType The type of problem being solved.
     */
    template <typename ProblemType> struct IpIterate
    {
        typedef std::shared_ptr<Nlp<ProblemType>> NlpSp;

        /**
         * @brief Constructs an IpIterate object.
         * @param nlp Shared pointer to the NLP problem.
         */
        IpIterate(const NlpSp &nlp);

    public:
        // Problem information
        Scalar objective_scale = 1.; ///< Scaling factor for the objective function.

        /**
         * @brief Initializes the IpIterate object.
         */
        void initialize();

        /**
         * @brief Resets all evaluated quantities, marking them as not computed.
         */
        void reset_evaluated_quantities();

        /**
         * @brief Returns the problem information.
         * @return Constant reference to the ProblemInfo object.
         */
        const ProblemInfo<ProblemType> &info() const { return info_; }

        /**
         * @brief Returns the NLP problem.
         * @return Constant reference to the NLP shared pointer.
         */
        const NlpSp &nlp() const { return nlp_; }

        /**
         * @brief Returns the primal variables x.
         * @return Constant reference to the primal x vector.
         */
        const VecRealView &primal_x() const { return primal_x_; }

        /**
         * @brief Returns the primal slack variables s.
         * @return Constant reference to the primal s vector.
         */
        const VecRealView &primal_s() const { return primal_s_; }

        /**
         * @brief Returns the dual variables for equality constraints.
         * @return Constant reference to the dual equality vector.
         */
        const VecRealView &dual_eq() const { return dual_eq_; }

        /**
         * @brief Returns the dual variables for lower bound constraints.
         * @return Constant reference to the dual lower bounds vector.
         */
        const VecRealView &dual_bounds_l() const { return dual_bounds_l_; }

        /**
         * @brief Returns the dual variables for upper bound constraints.
         * @return Constant reference to the dual upper bounds vector.
         */
        const VecRealView &dual_bounds_u() const { return dual_bounds_u_; }

        /**
         * @brief Returns the search direction for primal variables x.
         * @return Constant reference to the delta primal x vector.
         */
        const VecRealView &delta_primal_x() const { return delta_primal_x_; }

        /**
         * @brief Returns the search direction for primal slack variables s.
         * @return Constant reference to the delta primal s vector.
         */
        const VecRealView &delta_primal_s() const { return delta_primal_s_; }

        /**
         * @brief Returns the search direction for dual equality variables.
         * @return Constant reference to the delta dual equality vector.
         */
        const VecRealView &delta_dual_eq() const { return delta_dual_eq_; }

        /**
         * @brief Returns the search direction for dual lower bound variables.
         * @return Constant reference to the delta dual lower bounds vector.
         */
        const VecRealView &delta_dual_bounds_l() const { return delta_dual_bounds_l_; }

        /**
         * @brief Returns the search direction for dual upper bound variables.
         * @return Constant reference to the delta dual upper bounds vector.
         */
        const VecRealView &delta_dual_bounds_u() const { return delta_dual_bounds_u_; }

        /**
         * @brief Returns the lower bounds of the variables.
         * @return Constant reference to the lower bounds vector.
         */
        const VecRealView &lower_bounds() const { return lower_bounds_; };

        /**
         * @brief Returns the upper bounds of the variables.
         * @return Constant reference to the upper bounds vector.
         */
        const VecRealView &upper_bounds() const { return upper_bounds_; };

        /**
         * @brief Returns a boolean vector indicating which variables are lower bounded.
         * @return Constant reference to the lower bounded boolean vector.
         */
        const std::vector<bool> &lower_bounded() const { return lower_bounded_; };

        /**
         * @brief Returns a boolean vector indicating which variables are upper bounded.
         * @return Constant reference to the upper bounded boolean vector.
         */
        const std::vector<bool> &upper_bounded() const { return upper_bounded_; };

        /**
         * @brief Sets the primal variables x.
         * @param primal_x The new primal x vector.
         */
        template <typename Derived> void set_primal_x(const VecReal<Derived> &primal_x);

        /**
         * @brief Sets the primal slack variables s.
         * @param primal_s The new primal s vector.
         */
        template <typename Derived> void set_primal_s(const VecReal<Derived> &primal_s);

        /**
         * @brief Sets the dual variables for equality constraints.
         * @param dual_eq The new dual equality vector.
         */
        template <typename Derived> void set_dual_eq(const VecReal<Derived> &dual_eq);

        /**
         * @brief Sets the dual variables for lower bound constraints.
         * @param dual_bounds_l The new dual lower bounds vector.
         */
        template <typename Derived> void set_dual_bounds_l(const VecReal<Derived> &dual_bounds_l);

        /**
         * @brief Sets the dual variables for upper bound constraints.
         * @param dual_bounds_u The new dual upper bounds vector.
         */
        template <typename Derived> void set_dual_bounds_u(const VecReal<Derived> &dual_bounds_u);

        /**
         * @brief Sets the search direction for primal variables x.
         * @param delta_primal_x The new delta primal x vector.
         */
        template <typename Derived> void set_delta_primal_x(const VecReal<Derived> &delta_primal_x);

        /**
         * @brief Sets the search direction for primal slack variables s.
         * @param delta_primal_s The new delta primal s vector.
         */
        template <typename Derived> void set_delta_primal_s(const VecReal<Derived> &delta_primal_s);

        /**
         * @brief Sets the search direction for dual equality variables.
         * @param delta_dual_eq The new delta dual equality vector.
         */
        template <typename Derived> void set_delta_dual_eq(const VecReal<Derived> &delta_dual_eq);

        /**
         * @brief Sets the search direction for dual lower bound variables.
         * @param delta_dual_bounds_l The new delta dual lower bounds vector.
         */
        template <typename Derived>
        void set_delta_dual_bounds_l(const VecReal<Derived> &delta_dual_bounds_l);

        /**
         * @brief Sets the search direction for dual upper bound variables.
         * @param delta_dual_bounds_u The new delta dual upper bounds vector.
         */
        template <typename Derived>
        void set_delta_dual_bounds_u(const VecReal<Derived> &delta_dual_bounds_u);

        /**
         * @brief Computes and returns the objective function value.
         * @return The objective function value.
         */
        Scalar obj_value();

        /**
         * @brief Computes and returns the barrier function value.
         * @return The barrier function value.
         */
        Scalar barrier_value();

        /**
         * @brief Returns the gradient of the objective function with respect to x.
         * @return Constant reference to the objective gradient w.r.t. x.
         */
        const VecRealView &obj_gradient_x();

        /**
         * @brief Returns the gradient of the objective function with respect to s.
         * @return Constant reference to the objective gradient w.r.t. s.
         */
        const VecRealView &obj_gradient_s();

        /**
         * @brief Computes and returns the constraint violation vector.
         * @return Constant reference to the constraint violation vector.
         */
        const VecRealView &constr_viol();

        /**
         * @brief Computes and returns the inequality constraint violation vector.
         * @return Constant reference to the inequality constraint violation vector.
         * @note This requires specialization for the ProblemType.
         */
        const VecRealView constr_viol_ineq();

        /**
         * @brief Computes and returns the dual infeasibility with respect to x.
         * @return Constant reference to the dual infeasibility w.r.t. x.
         */
        const VecRealView &dual_infeasibility_x();

        /**
         * @brief Computes and returns the dual infeasibility with respect to s.
         * @return Constant reference to the dual infeasibility w.r.t. s.
         */
        const VecRealView &dual_infeasibility_s();

        /**
         * @brief Computes and returns the gradient of the barrier function.
         * @return Constant reference to the barrier function gradient.
         */
        const VecRealView &barrier_gradient();

        /**
         * @brief Computes and returns s - L when lower bounded, else 1.
         * @return Constant reference to the delta_lower vector.
         */
        const VecRealView &delta_lower();

        /**
         * @brief Computes and returns U - s when upper bounded, else 1.
         * @return Constant reference to the delta_upper vector.
         */
        const VecRealView &delta_upper();

        /**
         * @brief Computes and returns the complementarity for lower bounds.
         * @return Constant reference to the complementarity vector for lower bounds.
         */
        const VecRealView &complementarity_l();

        /**
         * @brief Computes and returns the complementarity for upper bounds.
         * @return Constant reference to the complementarity vector for upper bounds.
         */
        const VecRealView &complementarity_u();

        /**
         * @brief Returns the current value of mu (barrier parameter).
         * @return The current value of mu.
         */
        Scalar mu() const { return mu_; };

        /**
         * @brief Sets the value of mu (barrier parameter).
         * @param value The new value for mu.
         */
        void set_mu(Scalar value);

        /**
         * @brief Returns the primal inertia correction vector Dx_.
         * @return Constant reference to the Dx_ vector.
         */
        const VecRealView &Dx() const { return Dx_; }

        /**
         * @brief Returns the equality inertia correction vector De_.
         * @return Constant reference to the De_ vector.
         */
        const VecRealView &De() const { return De_; }

        /**
         * @brief Returns the flag indicating if De_ is zero.
         * @return The value of De_is_zero_.
         */
        bool De_is_zero() const { return De_is_zero_; }

        /**
         * @brief Sets the primal inertia correction vector Dx_.
         * @param Dx The new Dx_ vector.
         */
        template <typename Derived> void set_Dx(const VecReal<Derived> &Dx);

        /**
         * @brief Sets the equality inertia correction vector De_.
         * @param De The new De_ vector.
         */
        template <typename Derived> void set_De(const VecReal<Derived> &De);

        /**
         * @brief Sets the flag indicating if De_ is zero.
         * @param value The new value for De_is_zero_.
         */
        void set_De_is_zero(bool value);

        /**
         * @brief Computes the linear decrease in the objective function.
         * @return The linear decrease in the objective function.
         */
        Scalar linear_decrease_objective();

        /**
         * @brief Computes the linear decrease in the barrier function.
         * @return The linear decrease in the barrier function.
         */
        Scalar linear_decrease_barrier();

        /**
         * @brief Computes the error for a given mu value.
         * @param mu The mu value to compute the error for.
         * @return The computed error.
         */
        Scalar e_mu(Scalar mu);

        Scalar tau() const { return std::max(0.99, 1.-mu()); }

        /**
         * @brief Computes the maximum step size.
         * @param tau The fraction-to-the-boundary parameter.
         * @return A pair of Scalars representing the maximum primal and dual step sizes.
         */
        std::pair<Scalar, Scalar> maximum_step_size(const Scalar tau);

        /**
         * @brief Returns the number of bounds in the problem.
         * @return The number of bounds.
         */
        Index number_of_bounds() const { return number_of_bounds_; }

        /**
         * @brief Returns a reference to the Hessian of the problem.
         * @return Reference to the Hessian.
         */
        Hessian<ProblemType> &hessian();

        /**
         * @brief Returns a reference to a zero Hessian.
         * @return Reference to a zero Hessian.
         */
        Hessian<ProblemType> &zero_hessian();

        /**
         * @brief Returns a reference to the Jacobian of the problem.
         * @return Reference to the Jacobian.
         */
        Jacobian<ProblemType> &jacobian();

    private:
        const ProblemInfo<ProblemType> info_; ///< Information about the NLP.
        NlpSp nlp_;
        // Iteration point
        VecRealAllocated primal_x_;      ///< Primal variables of the NLP.
        VecRealAllocated primal_s_;      ///< Primal variables of the NLP.
        VecRealAllocated dual_eq_;       ///< Dual variables of the equality constraints.
        VecRealAllocated dual_bounds_l_; ///< Dual variables of the bound constraints.
        VecRealAllocated dual_bounds_u_; ///< Dual variables of the bound constraints.
        // Searh direction
        VecRealAllocated delta_primal_x_; ///< Primal search direction.
        VecRealAllocated delta_primal_s_; ///< Primal search direction.
        VecRealAllocated delta_dual_eq_;  ///< Dual search direction of the equality constraints.
        VecRealAllocated
            delta_dual_bounds_l_; ///< Dual search direction of the lower bound constraints.
        VecRealAllocated
            delta_dual_bounds_u_; ///< Dual search direction of the upper bound constraints.
        // Evaluated quantities
        Scalar obj_value_;                      ///< Objective value of the NLP.
        Scalar barrier_value_;                  ///< Value of the barrier.
        VecRealAllocated obj_gradient_x_;       ///< Gradient of the objective function.
        VecRealAllocated obj_gradient_s_;       ///< Gradient of the objective function.
        VecRealAllocated constr_viol_;          ///< Residual of the equality constraints.
        VecRealAllocated dual_infeasibility_x_; ///< Residual of the equality constraints.
        VecRealAllocated dual_infeasibility_s_; ///< Residual of the equality constraints.
        // Derived quantities
        Scalar mu_;                         ///< Barrier value of the NLP.
        VecRealAllocated barrier_gradient_; ///< Gradient of the barrier function.
        VecRealAllocated delta_lower_;
        VecRealAllocated delta_upper_;
        Scalar linear_decrease_objective_;   ///< Linear decrease of the objective function.
        Scalar linear_decrease_barrier_;     ///< Linear decrease of the barrier function.
        VecRealAllocated complementarity_l_; ///< Complementarity of the NLP.
        VecRealAllocated complementarity_u_; ///< Complementarity of the NLP.
        Jacobian<ProblemType> jacobian_;     ///< Jacobian of the NLP.
        Hessian<ProblemType> hessian_;       ///< Hessian of the NLP.
        // Problem information
        VecRealAllocated lower_bounds_; ///< Lower bounds of the variables.
        VecRealAllocated upper_bounds_; ///< Upper bounds of the variables.
        std::vector<bool>
            lower_bounded_; ///< Boolean vector indicating if the variables are lower bounded.
        std::vector<bool>
            upper_bounded_; ///< Boolean vector indicating if the variables are upper bounded.
        std::vector<bool> single_lower_bounded_; ///< Boolean vector indicating if the variables are
                                                 ///< lower bounded.
        std::vector<bool> single_upper_bounded_; ///< Boolean vector indicating if the variables are
                                                 ///< upper bounded.
        // Vector quantities related to solving for the search direction
        VecRealAllocated Dx_; ///< primal inertia correction vector
        VecRealAllocated De_; ///< equality inertia correction vector
        bool De_is_zero_ = false;

        // Flags to track which quantities have already been evaluated
        bool obj_value_evaluated_ = false;      ///< Flag for objective value evaluation
        bool obj_gradient_evaluated_ = false;   ///< Flag for objective gradient (x) evaluation
        bool obj_gradient_evaluated_s_ = false; ///< Flag for objective gradient (s) evaluation
        bool constr_viol_evaluated_ = false;    ///< Flag for constraint violation evaluation
        bool dual_infeasibility_x_evaluated_ =
            false; ///< Flag for dual infeasibility (x) evaluation
        bool dual_infeasibility_s_evaluated_ =
            false;                                ///< Flag for dual infeasibility (s) evaluation
        bool jacobian_evaluated_ = false;         ///< Flag for Jacobian evaluation
        bool hessian_evaluated_ = false;          ///< Flag for Hessian evaluation
        bool barrier_value_evaluated_ = false;    ///< Flag for barrier value evaluation
        bool barrier_gradient_evaluated_ = false; ///< Flag for barrier gradient evaluation
        bool delta_lower_evaluated_ = false;      ///< Flag for delta_lower evaluation
        bool delta_upper_evaluated_ = false;      ///< Flag for delta_upper evaluation
        bool linear_decrease_objective_evaluated_ =
            false; ///< Flag for linear decrease in objective evaluation
        bool linear_decrease_barrier_evaluated_ =
            false; ///< Flag for linear decrease in barrier evaluation
        bool complementarity_l_evaluated_ =
            false; ///< Flag for lower bound complementarity evaluation
        bool complementarity_u_evaluated_ =
            false; ///< Flag for upper bound complementarity evaluation

        Scalar kappa_d_ = 0.;      ///< Fraction-to-the-boundary parameter
        Index number_of_bounds_ = 0; ///< Total number of bounds in the problem
        Scalar smax_ = 100.;         ///< Maximum allowed slack variable value
    };
} // namespace fatrop
// implementation
template <typename ProblemType>
template <typename Derived>
void fatrop::IpIterate<ProblemType>::set_primal_x(const VecReal<Derived> &primal_x)
{
    reset_evaluated_quantities();
    primal_x_ = primal_x;
}
template <typename ProblemType>
template <typename Derived>
void fatrop::IpIterate<ProblemType>::set_primal_s(const VecReal<Derived> &primal_s)
{
    reset_evaluated_quantities();
    primal_s_ = primal_s;
}
template <typename ProblemType>
template <typename Derived>
void fatrop::IpIterate<ProblemType>::set_dual_eq(const VecReal<Derived> &dual_eq)
{
    dual_infeasibility_s_evaluated_ = false;
    dual_infeasibility_x_evaluated_ = false;
    dual_eq_ = dual_eq;
}
template <typename ProblemType>
template <typename Derived>
void fatrop::IpIterate<ProblemType>::set_dual_bounds_l(const VecReal<Derived> &dual_bounds_l)
{
    complementarity_l_evaluated_ = false;
    dual_infeasibility_s_evaluated_ = false;
    dual_bounds_l_ = dual_bounds_l;
}
template <typename ProblemType>
template <typename Derived>
void fatrop::IpIterate<ProblemType>::set_dual_bounds_u(const VecReal<Derived> &dual_bounds_u)
{
    complementarity_u_evaluated_ = false;
    dual_infeasibility_s_evaluated_ = false;
    dual_bounds_u_ = dual_bounds_u;
}
template <typename ProblemType>
template <typename Derived>
void fatrop::IpIterate<ProblemType>::set_delta_primal_x(const VecReal<Derived> &delta_primal_x)
{
    linear_decrease_objective_evaluated_ = false;
    delta_primal_x_ = delta_primal_x;
}
template <typename ProblemType>
template <typename Derived>
void fatrop::IpIterate<ProblemType>::set_delta_primal_s(const VecReal<Derived> &delta_primal_s)
{
    linear_decrease_barrier_evaluated_ = false;
    linear_decrease_objective_evaluated_ = false;
    delta_primal_s_ = delta_primal_s;
}
template <typename ProblemType>
template <typename Derived>
void fatrop::IpIterate<ProblemType>::set_delta_dual_eq(const VecReal<Derived> &delta_dual_eq)
{
    delta_dual_eq_ = delta_dual_eq;
}
template <typename ProblemType>
template <typename Derived>
void fatrop::IpIterate<ProblemType>::set_delta_dual_bounds_l(
    const VecReal<Derived> &delta_dual_bounds_l)
{
    delta_dual_bounds_l_ = delta_dual_bounds_l;
}
template <typename ProblemType>
template <typename Derived>
void fatrop::IpIterate<ProblemType>::set_delta_dual_bounds_u(
    const VecReal<Derived> &delta_dual_bounds_u)
{
    delta_dual_bounds_u_ = delta_dual_bounds_u;
}

template <typename ProblemType>
template <typename Derived>
void fatrop::IpIterate<ProblemType>::set_Dx(const VecReal<Derived> &Dx)
{
    Dx_ = Dx;
}

template <typename ProblemType>
template <typename Derived>
void fatrop::IpIterate<ProblemType>::set_De(const VecReal<Derived> &De)
{
    De_ = De;
}

template <typename ProblemType>
void fatrop::IpIterate<ProblemType>::set_De_is_zero(bool value)
{
    De_is_zero_ = value;
}

#endif //__fatrop_algorithm_ip_iterate_hpp__
