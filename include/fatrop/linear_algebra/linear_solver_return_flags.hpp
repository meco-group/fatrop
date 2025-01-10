#ifndef __fatrop_linear_algebra_linear_solver_return_flags_hpp__
#define __fatrop_linear_algebra_linear_solver_return_flags_hpp__
/**
 * @enum LinsolReturnFlag
 * @brief Enumeration of possible return flags for the linear solver.
 */
namespace fatrop
{
    enum LinsolReturnFlag
    {
        SUCCESS = 0,     ///< The solver successfully found a solution.
        INDEFINITE = 1,  ///< The reduced Hessian is indefinite; no descent direction found.
        NOFULL_RANK = 2, ///< The Jacobian is (numerically) not full row rank.
        ITREF_MAX_ITER =
            3, ///<  Iterative refinement did not converge within the maximum number of iterations.
        ITREF_INCREASE = 4, ///< Iterative refinement did not converge due to increased residual.
        UNKNOWN = 5         ///< An unknown error occurred during the solving process.
    };
} // namespace fatrop
#endif //__fatrop_linear_algebra_linear_solver_return_flags_hpp__