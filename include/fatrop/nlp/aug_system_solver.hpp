//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_nlp_aug_system_solver_hpp__
#define __fatrop_nlp_aug_system_solver_hpp__
#include "fatrop/nlp/fwd.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"
#include "fatrop/linear_algebra/fwd.hpp"

namespace fatrop
{
    template<typename ProblemType>
    class AugSystemSolver
    {
        /**
         * @brief Solves the augmented system without path equality constraint regularization.
         * @param info Problem information.
         * @param jacobian Jacobian of the constraints.
         * @param hessian Hessian of the Lagrangian.
         * @param D_x Diagonal regularization for primal variables.
         * @param D_s Diagonal regularization for slack variables.
         * @param f Gradient of the objective function.
         * @param g Constraint residuals.
         * @param x [out] Solution vector for primal variables.
         * @param eq_mult [out] Solution vector for equality constraint multipliers.
         * @return Status flag indicating the outcome of the solve operation.
         */
        LinsolReturnFlag solve(const ProblemInfo<ProblemType> &info, Jacobian<ProblemType> &jacobian,
                               Hessian<ProblemType> &hessian, const VecRealView &D_x,
                               const VecRealView &D_s, const VecRealView &f, const VecRealView &g,
                               VecRealView &x, VecRealView &eq_mult);

        /**
         * @brief Solves the augmented system with path equality constraint regularization.
         * @param info Problem information.
         * @param jacobian Jacobian of the constraints.
         * @param hessian Hessian of the Lagrangian.
         * @param D_x Diagonal regularization for primal variables.
         * @param D_eq Diagonal regularization for equality constraints.
         * @param D_s Diagonal regularization for slack variables.
         * @param f Gradient of the objective function.
         * @param g Constraint residuals.
         * @param x [out] Solution vector for primal variables.
         * @param eq_mult [out] Solution vector for equality constraint multipliers.
         * @return Status flag indicating the outcome of the solve operation.
         */
        LinsolReturnFlag solve(const ProblemInfo<ProblemType> &info, Jacobian<ProblemType> &jacobian,
                               Hessian<ProblemType> &hessian, const VecRealView &D_x,
                               const VecRealView &D_eq, const VecRealView &D_s,
                               const VecRealView &f, const VecRealView &g, VecRealView &x,
                               VecRealView &eq_mult);

        /**
         * @brief Solves the system for a new right-hand side without path equality constraint regularization.
         * @param info Problem information.
         * @param jacobian Jacobian of the constraints.
         * @param hessian Hessian of the Lagrangian.
         * @param D_s Diagonal regularization for slack variables.
         * @param f Gradient of the objective function.
         * @param g Constraint residuals.
         * @param x [out] Solution vector for primal variables.
         * @param eq_mult [out] Solution vector for equality constraint multipliers.
         * @return Status flag indicating the outcome of the solve operation.
         */
        LinsolReturnFlag solve_rhs(const ProblemInfo<ProblemType> &info,
                                   const Jacobian<ProblemType> &jacobian,
                                   const Hessian<ProblemType> &hessian, const VecRealView &D_s,
                                   const VecRealView &f, const VecRealView &g, VecRealView &x,
                                   VecRealView &eq_mult);

        /**
         * @brief Solves the system for a new right-hand side with path equality constraint regularization.
         * @param info Problem information.
         * @param jacobian Jacobian of the constraints.
         * @param hessian Hessian of the Lagrangian.
         * @param D_eq Diagonal regularization for equality constraints.
         * @param D_s Diagonal regularization for slack variables.
         * @param f Gradient of the objective function.
         * @param g Constraint residuals.
         * @param x [out] Solution vector for primal variables.
         * @param eq_mult [out] Solution vector for equality constraint multipliers.
         * @return Status flag indicating the outcome of the solve operation.
         */
        LinsolReturnFlag solve_rhs(const ProblemInfo<ProblemType> &info,
                                   const Jacobian<ProblemType> &jacobian,
                                   const Hessian<ProblemType> &hessian, const VecRealView &D_eq,
                                   const VecRealView &D_s, const VecRealView &f,
                                   const VecRealView &g, VecRealView &x, VecRealView &eq_mult);
    };

} // namespace fatrop

#endif //__fatrop_nlp_problem_info__
