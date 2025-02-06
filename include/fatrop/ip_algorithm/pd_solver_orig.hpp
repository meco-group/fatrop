//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_pd_solver_orig__
#define __fatrop_ip_algorithm_pd_solver_orig__

#include "fatrop/linear_algebra/linear_solver.hpp"
#include "fatrop/ip_algorithm/pd_system_orig.hpp"
#include "fatrop/nlp/fwd.hpp"

namespace fatrop
{
    /**
     * @brief Primal-dual solver for the original problem formulation.
     * 
     * This class is intended to implement a linear solver for the primal-dual system
     * arising in interior point methods for the original problem formulation.
     * 
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template<typename ProblemType>
    class PdSolverOrig;
    // : public LinearSolver<PdSolverOrig<ProblemType>, PdSystemType<ProblemType>>
    // {
    //     /**
    //      * @brief Solve the linear system once.
    //      * @param ls The primal-dual system to solve.
    //      * @param x The solution vector.
    //      * @return LinsolReturnFlag Indicating the success or failure of the solve.
    //      */
    //     LinsolReturnFlag solve_once_impl(PdSystemType<ProblemType> &ls, VecRealView &x);

    //     /**
    //      * @brief Solve the system for a new right-hand side.
    //      * @param ls The primal-dual system to solve.
    //      * @param x The solution vector.
    //      */
    //     void solve_rhs_impl(PdSystemType<ProblemType> &ls, VecRealView &x);
    // };
} // namespace fatrop

#endif //__fatrop_ip_algorithm_pd_solver_orig__
