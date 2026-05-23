//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_linear_algebra_linear_solver_hpp__
#define __fatrop_linear_algebra_linear_solver_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"
#include "fatrop/linear_algebra/vector.hpp"

namespace fatrop
{
    /**
     * @brief A template class for linear solvers.
     *
     * @tparam Derived The derived class type (CRTP pattern).
     * @tparam LsType The linear system type.
     */
    template <typename Derived, typename LsType> class LinearSolver
    {
    public:
        typedef LinearSystem<LsType> LinearSystemType;
        /**
         * @brief Construct a new Linear Solver object.
         *
         * @param m The dimension of the linear system.
         */
        LinearSolver(const Index m) : m_(m), x(m), x_(m_), residual_(m_), tmp_(m_) {}
        /**
         * @brief Get a reference to the derived class.
         *
         * @return Derived& Reference to the derived class.
         */
        Derived &derived() { return *static_cast<Derived *>(this); }

        /**
         * @brief Solve the linear system without iterative refinement.
         *
         * @param ls The linear system to solve.
         * @param x_in The input/output vector for the solution.
         * @return LinsolReturnFlag Indicates the success or failure of the solve operation.
         */
        LinsolReturnFlag solve_once(LinearSystem<LsType> &ls, VecRealView &x_in);

        /**
         * @brief Solve the right-hand side of the linear system.
         *
         * @param ls The linear system.
         * @param x_in The input/output vector for the solution.
         */
        void solve_rhs(LinearSystem<LsType> &ls, VecRealView &x_in);

        /**
         * @brief Solve the linear system with iterative refinement.
         *
         * @param ls The linear system to solve.
         * @return LinsolReturnFlag Indicates the success or failure of the solve operation.
         */
        LinsolReturnFlag solve_in_place(LinearSystem<LsType> &ls);

        /**
         * @brief Solve the right-hand side of the linear system with iterative refinement.
         *
         * This method assumes that solve_in_place had been called before, such that the
         * factorization can be reused.
         *
         * @param ls The linear system.
         * @return LinsolReturnFlag Indicates the success or failure of the solve operation.
         */
        LinsolReturnFlag solve_in_place_rhs(LinearSystem<LsType> &ls);


        // Iterative-refinement tuning + monitoring.
        void set_min_it_ref(const Index &value) { min_it_ref = value; }
        void set_max_it_ref(const Index &value) { max_it_ref = value; }
        void set_iref_tol(const Scalar &value) { tol_ = value; }
        /// Last (final) residual norm = ||A x + b||_inf / max(||b||, 1)
        /// produced by `apply_iterative_refinement`.
        Scalar last_iref_residual() const { return last_residual_; }
        /// Number of iterative-refinement iterations actually performed by
        /// the last call to `apply_iterative_refinement`.
        Index last_iref_iters() const { return last_iref_iters_; }
        /// Largest iterative-refinement residual seen since construction
        /// (or since @c reset_iref_stats()).
        Scalar worst_iref_residual() const { return worst_residual_; }
        void reset_iref_stats() { worst_residual_ = 0.0; }

    protected:
        /**
         * @brief Apply iterative refinement to improve the solution accuracy.
         *
         * @param ls The linear system.
         * @return LinsolReturnFlag Indicates the success or failure of the refinement process.
         */
        LinsolReturnFlag apply_iterative_refinement(LinearSystem<LsType> &ls);

        const Index m_;             ///< Dimension of the linear system.
        Index min_it_ref = 0;       ///< Minimum number of iterative refinement steps.
        Index max_it_ref = 5;      ///< Maximum number of iterative refinement steps.
        Scalar tol_ = 1e-8;         ///< Tolerance for iterative refinement convergence.
        VecRealAllocated x;         ///< Allocated vector for solution.
        VecRealAllocated x_;        ///< Allocated vector for intermediate solution.
        VecRealAllocated residual_; ///< Allocated vector for residual.
        VecRealAllocated tmp_;      ///< Allocated vector for temporary calculations.
        // Monitoring stats (updated by apply_iterative_refinement).
        Scalar last_residual_ = 0.0;
        Index last_iref_iters_ = 0;
        Scalar worst_residual_ = 0.0;
    };
} // namespace fatrop

#endif // __fatrop_linear_algebra_linear_solver_hpp__
