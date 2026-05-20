//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_nlp_dims__
#define __fatrop_nlp_dims__

#include "fatrop/context/context.hpp"

namespace fatrop
{
    struct NlpDims
    {
        // Backwards-compatible constructor: the tangent (search-direction) dimension
        // defaults to the primal (manifold) dimension, recovering the standard
        // Euclidean behavior.
        NlpDims(const Index number_of_variables, const Index number_of_eq_constraints,
                const Index number_of_ineq_constraints)
            : number_of_variables(number_of_variables),
              number_of_tangent_variables(number_of_variables),
              number_of_eq_constraints(number_of_eq_constraints),
              number_of_ineq_constraints(number_of_ineq_constraints) {};
        // Lie-group / manifold variant: primal and tangent dimensions differ.
        NlpDims(const Index number_of_variables, const Index number_of_tangent_variables,
                const Index number_of_eq_constraints, const Index number_of_ineq_constraints)
            : number_of_variables(number_of_variables),
              number_of_tangent_variables(number_of_tangent_variables),
              number_of_eq_constraints(number_of_eq_constraints),
              number_of_ineq_constraints(number_of_ineq_constraints) {};
        NlpDims() {};
        Index number_of_variables;         ///< Size of the primal (manifold) variable vector x.
        Index number_of_tangent_variables; ///< Size of the search direction / tangent vector
                                           ///< delta_x. The Newton system, gradients and the
                                           ///< Jacobian/Hessian column/row blocks live in this
                                           ///< space. Defaults to number_of_variables for
                                           ///< standard Euclidean problems.
        Index number_of_eq_constraints;    ///< Number of equality constraints in the NLP.
        Index number_of_ineq_constraints;  ///< Number of inequality constraints in the NLP.
        /**
         * todo: rename to number of slack variables
         */
    };
    template <typename ProblemType> struct ProblemDims
    {
    };

} // namespace fatrop

#endif //__fatrop_nlp_dims__
