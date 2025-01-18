//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_algorithm_ip_iterate_hpp__
#define __fatrop_algorithm_ip_iterate_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/nlp/hessian.hpp"
#include "fatrop/nlp/jacobian.hpp"
#include "fatrop/nlp/nlp.hpp"
#include <memory>

namespace fatrop
{
    template <typename ProblemType> struct IpIterate
    {
        typedef std::shared_ptr<Nlp<ProblemType>> NlpSp;
        IpIterate(const NlpSp &nlp);

    public:
        // Problem information
        ProblemInfo<ProblemType> info; ///< Information about the NLP.
        // Iteration point
        VecRealAllocated primal;        ///< Primal variables of the NLP.
        VecRealAllocated primal_s;      ///< Primal variables of the NLP.
        VecRealAllocated dual_eq;       ///< Dual variables of the equality constraints.
        VecRealAllocated dual_bounds_L; ///< Dual variables of the bound constraints.
        VecRealAllocated dual_bounds_Z; ///< Dual variables of the bound constraints.
        // Searh direction
        VecRealAllocated delta_primal;   ///< Primal search direction.
        VecRealAllocated delta_primal_s; ///< Primal search direction.
        VecRealAllocated delta_dual_eq;  ///< Dual search direction of the equality constraints.
        VecRealAllocated
            delta_dual_bounds_L; ///< Dual search direction of the lower bound constraints.
        VecRealAllocated
            delta_dual_bounds_Z; ///< Dual search direction of the upper bound constraints.
        // Evaluated quantities
        Scalar obj_value;               ///< Objective value of the NLP.
        VecRealAllocated obj_gradient;  ///< Gradient of the objective function.
        VecRealAllocated constr_viol;   ///< Residual of the equality constraints.
        VecRealAllocated dual_eq_feas;  ///< Residual of the equality constraints.
        // Derived quantities
        Scalar barrier_value;              ///< Barrier value of the NLP.
        VecRealAllocated barrier_gradient; ///< Gradient of the barrier function.
        Scalar linear_decrease;            ///< Linear decrease of the barrier function.
        Scalar objective_scale = 1.;       ///< Scaling factor for the objective function.
        void reset_evaluated_quantities()
        {
            obj_value_evaluated_ = false;
            obj_gradient_evaluated_ = false;
            constr_viol_evaluated_ = false;
            dual_eq_feas_evaluated_ = false;
            jacobian_evaluated_ = false;
            hessian_evaluated_ = false;
            barrier_value_evaluated_ = false;
            barrier_gradient_evaluated_ = false;
            linear_decrease_evaluated_ = false;
        }

        Hessian<ProblemType> &hessian();
        Jacobian<ProblemType> &jacobian();
        // VecRealView &constr_viol();

    private:
        NlpSp nlp_;
        Jacobian<ProblemType> jacobian_; ///< Jacobian of the NLP.
        Hessian<ProblemType> hessian_;   ///< Hessian of the NLP.
        // booleans to keep track of which quantities have already been evaluated
        bool obj_value_evaluated_ = false;
        bool obj_gradient_evaluated_ = false;
        bool constr_viol_evaluated_ = false;
        bool dual_eq_feas_evaluated_ = false;
        bool jacobian_evaluated_ = false;
        bool hessian_evaluated_ = false;
        bool barrier_value_evaluated_ = false;
        bool barrier_gradient_evaluated_ = false;
        bool linear_decrease_evaluated_ = false;
    };
} // namespace fatrop

#endif //__fatrop_algorithm_ip_iterate_hpp__