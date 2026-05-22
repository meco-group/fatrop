//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_dense_type_hpp__
#define __fatrop_dense_type_hpp__

/**
 * @file type.hpp
 * @brief Defines the DenseType structure for dense NLP problems.
 *
 * This file contains the declaration of the DenseType structure,
 * which is used as a base for problem class-specific specializations
 * (Jacobian, Hessian, AugSystemSolver, ...) for dense problems.
 *
 * A dense problem has the form
 *     min  f(x)
 *     s.t. g_eq(x)   = 0
 *          g_ineq(x) in [L, U]
 * with a single dense primal vector x of size nx, no time/staging
 * structure and no separate control variables. The DenseType linear
 * solver is morally equivalent to running the OcpType solver with
 * K = 1, nu = 0, but the code is specialized and simplified.
 */

namespace fatrop
{
    struct DenseType
    {
    };
} // namespace fatrop

#endif // __fatrop_dense_type_hpp__
