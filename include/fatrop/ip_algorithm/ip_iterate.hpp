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
    template <typename ProblemType> struct IpIterate
    {
        typedef std::shared_ptr<Nlp<ProblemType>> NlpSp;
        IpIterate(const NlpSp &nlp);

    public:
        // Problem information
        Scalar objective_scale = 1.;   ///< Scaling factor for the objective function.
        void initialize();
        void reset_evaluated_quantities();

        const ProblemInfo<ProblemType> &info() const { return info_; }
        const NlpSp &nlp() const { return nlp_; }

        const VecRealView &primal_x() const { return primal_x_; }
        const VecRealView &primal_s() const { return primal_s_; }
        const VecRealView &dual_eq() const { return dual_eq_; }
        const VecRealView &dual_bounds_l() const { return dual_bounds_l_; }
        const VecRealView &dual_bounds_u() const { return dual_bounds_u_; }
        const VecRealView &delta_primal_x() const { return delta_primal_x_; }
        const VecRealView &delta_primal_s() const { return delta_primal_s_; }
        const VecRealView &delta_dual_eq() const { return delta_dual_eq_; }
        const VecRealView &delta_dual_bounds_l() const { return delta_dual_bounds_l_; }
        const VecRealView &delta_dual_bounds_u() const { return delta_dual_bounds_u_; }

        template <typename Derived> void set_primal_x(const VecReal<Derived> &primal_x);
        template <typename Derived> void set_primal_s(const VecReal<Derived> &primal_s);
        template <typename Derived> void set_dual_eq(const VecReal<Derived> &dual_eq);
        template <typename Derived> void set_dual_bounds_l(const VecReal<Derived> &dual_bounds_l);
        template <typename Derived> void set_dual_bounds_u(const VecReal<Derived> &dual_bounds_u);
        template <typename Derived> void set_delta_primal_x(const VecReal<Derived> &delta_primal_x);
        template <typename Derived> void set_delta_primal_s(const VecReal<Derived> &delta_primal_s);
        template <typename Derived> void set_delta_dual_eq(const VecReal<Derived> &delta_dual_eq);
        template <typename Derived>
        void set_delta_dual_bounds_l(const VecReal<Derived> &delta_dual_bounds_l);
        template <typename Derived>
        void set_delta_dual_bounds_u(const VecReal<Derived> &delta_dual_bounds_u);

        Scalar obj_value();
        Scalar barrier_value();
        const VecRealView &obj_gradient_x();
        const VecRealView &obj_gradient_s();
        const VecRealView &constr_viol();
        const VecRealView &dual_infeasibility_x();
        const VecRealView &dual_infeasibility_s();
        const VecRealView &barrier_gradient();
        // // s - L when lower bounded else 1.
        const VecRealView &delta_lower();
        // // U - s when upper bounded else 1.
        const VecRealView &delta_upper();
        const VecRealView &complementarity_l();
        const VecRealView &complementarity_u();
        Scalar mu() { return mu_; };
        void set_mu(Scalar value);
        Scalar linear_decrease_objective();
        Scalar linear_decrease_barrier();
        std::pair<Scalar, Scalar> maximum_step_size(const Scalar tau);

        Index number_of_bounds() const { return number_of_bounds_; }

        Hessian<ProblemType> &hessian();
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
        // booleans to keep track of which quantities have already been evaluated
        bool obj_value_evaluated_ = false;
        bool obj_gradient_evaluated_ = false;
        bool obj_gradient_evaluated_s_ = false;
        bool constr_viol_evaluated_ = false;
        bool dual_infeasibility_x_evaluated_ = false;
        bool dual_infeasibility_s_evaluated_ = false;
        bool jacobian_evaluated_ = false;
        bool hessian_evaluated_ = false;
        bool barrier_value_evaluated_ = false;
        bool barrier_gradient_evaluated_ = false;
        bool delta_lower_evaluated_ = false;
        bool delta_upper_evaluated_ = false;
        bool linear_decrease_objective_evaluated_ = false;
        bool linear_decrease_barrier_evaluated_ = false;
        bool complementarity_l_evaluated_ = false;
        bool complementarity_u_evaluated_ = false;
        Scalar kappa_d_ = 1e-5;
        Index number_of_bounds_ = 0;
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

#endif //__fatrop_algorithm_ip_iterate_hpp__