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
        NlpDims(const Index number_of_variables, const Index number_of_eq_constraints,
                const Index number_of_ineq_constraints)
            : number_of_variables(number_of_variables),
              number_of_eq_constraints(number_of_eq_constraints),
              number_of_ineq_constraints(number_of_ineq_constraints) {};
        const Index number_of_variables;        ///< Number of variables in the NLP.
        const Index number_of_eq_constraints;   ///< Number of equality constraints in the NLP.
        const Index number_of_ineq_constraints; ///< Number of inequality constraints in the NLP.
    };
    template <typename ProblemType> struct ProblemDims
    {
    };

} // namespace fatrop

#endif //__fatrop_nlp_dims__
