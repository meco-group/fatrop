//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_iterate_hpp__
#define __fatrop_ip_algorithm_ip_iterate_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/nlp/hessian.hpp"
#include "fatrop/nlp/jacobian.hpp"
#include "fatrop/nlp/fwd.hpp"
#include <memory>
#include <vector>
#include <utility>

namespace fatrop
{
    template <typename ProblemType> struct IpIterate
    {
        typedef std::shared_ptr<Nlp<ProblemType>> NlpSp;
        IpIterate(const NlpSp &nlp);

    public:
        // Problem information
        ProblemInfo<ProblemType> info; ///< Information about the NLP.
        Scalar objective_scale = 1.;   ///< Scaling factor for the objective function.
        void initialize();
        void reset_evaluated_quantities();

        VecRealView &primal_x() { return primal_x_; }
        const VecRealView &primal_x() const { return primal_x_; }
        VecRealView &primal_s() { return primal_s_; }
        const VecRealView &primal_s() const { return primal_s_; }
        VecRealView &dual_eq() { return dual_eq_; }
        const VecRealView &dual_eq() const { return dual_eq_; }
        VecRealView &dual_bounds_l() { return dual_bounds_l_; }
        const VecRealView &dual_bounds_l() const { return dual_bounds_l_; }
        VecRealView &dual_bounds_u() { return dual_bounds_u_; }
        const VecRealView &dual_bounds_u() const { return dual_bounds_u_; }
        VecRealView &delta_primal_x() { return delta_primal_x_; }
        const VecRealView &delta_primal_x() const { return delta_primal_x_; }
        VecRealView &delta_primal_s() { return delta_primal_s_; }
        const VecRealView &delta_primal_s() const { return delta_primal_s_; }
        VecRealView &delta_dual_eq() { return delta_dual_eq_; }
        const VecRealView &delta_dual_eq() const { return delta_dual_eq_; }
        VecRealView &delta_dual_bounds_l() { return delta_dual_bounds_l_; }
        const VecRealView &delta_dual_bounds_l() const { return delta_dual_bounds_l_; }
        VecRealView &delta_dual_bounds_u() { return delta_dual_bounds_u_; }
        const VecRealView &delta_dual_bounds_u() const { return delta_dual_bounds_u_; }

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

#endif //__fatrop_algorithm_ip_iterate_hpp__