//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_nlp_orig_hpp__
#define __fatrop_ip_algorithm_ip_nlp_orig_hpp__
#include "fatrop/common/fwd.hpp"
#include "fatrop/context/context.hpp"
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/nlp/fwd.hpp"
#include "fatrop/nlp/nlp.hpp"
#include <memory>
namespace fatrop
{
    /**
     * @brief Original NLP formulation for interior point methods.
     *
     * This class wraps the original NLP problem and provides the interface
     * required for the interior point algorithm.
     * Decorator class that adds some functionality to make the problem easier digestible by IP
     * algorithm.
     *
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename ProblemType> class IpNlpOrig : public Nlp<ProblemType>
    {
        typedef std::shared_ptr<Nlp<ProblemType>> NlpSp;

    public:
        /**
         * @brief Construct a new IpNlpOrig object.
         *
         * @param nlp Shared pointer to the original NLP problem.
         */
        IpNlpOrig(const NlpSp &nlp);

        /**
         * @brief Get the dimensions of the NLP problem.
         * @return const NlpDims& The dimensions of the NLP problem.
         */
        virtual const NlpDims &nlp_dims() const override;

        /**
         * @brief Get the dimensions of the specific problem type.
         * @return const ProblemDims<ProblemType>& The dimensions of the problem.
         */
        virtual const ProblemDims<ProblemType> &problem_dims() const override;

        /**
         * @brief Evaluate the Hessian of the Lagrangian.
         */
        virtual Index eval_lag_hess(const ProblemInfo<ProblemType> &info,
                                    const Scalar objective_scale, const VecRealView &primal_x,
                                    const VecRealView &primal_s, const VecRealView &lam,
                                    Hessian<ProblemType> &hess) override;

        /**
         * @brief Evaluate the Jacobian of the constraints.
         */
        virtual Index eval_constr_jac(const ProblemInfo<ProblemType> &info,
                                      const VecRealView &primal_x, const VecRealView &primal_s,
                                      Jacobian<ProblemType> &jac) override;

        /**
         * @brief Evaluate the constraint violation.
         */
        virtual Index eval_constraint_violation(const ProblemInfo<ProblemType> &info,
                                                const VecRealView &primal_x,
                                                const VecRealView &primal_s,
                                                VecRealView &res) override;

        /**
         * @brief Evaluate the gradient of the objective function.
         */
        virtual Index eval_objective_gradient(const ProblemInfo<ProblemType> &info,
                                              const Scalar objective_scale,
                                              const VecRealView &primal_x,
                                              const VecRealView &primal_s, VecRealView &grad_x,
                                              VecRealView &grad_s) override;

        /**
         * @brief Evaluate the objective function.
         */
        virtual Index eval_objective(const ProblemInfo<ProblemType> &info,
                                     const Scalar objective_scale, const VecRealView &primal_x,
                                     const VecRealView &primal_s, Scalar &res) override;

        /**
         * @brief Get the bounds of the problem.
         */
        virtual Index get_bounds(const ProblemInfo<ProblemType> &info, VecRealView &lower_bounds,
                                 VecRealView &upper_bounds) override;

        virtual Index get_initial_primal(const ProblemInfo<ProblemType> &info,
                                         VecRealView &primal_x) override;

        virtual void get_primal_damping(const ProblemInfo<ProblemType> &info,
                                        VecRealView &damping) override;

        virtual void apply_jacobian_s_transpose(const ProblemInfo<ProblemType> &info,
                                                const VecRealView &multipliers, const Scalar alpha,
                                                const VecRealView &y, VecRealView &out) override;
        virtual void callback(const IpData<ProblemType> &ip_data) override
        {
            nlp_->callback(ip_data);
        };

    private:
        NlpSp nlp_;                              ///< Shared pointer to the original NLP problem
        VecRealAllocated modified_bounds_lower_; ///< Modified lower bounds
        VecRealAllocated modified_bounds_upper_; ///< Modified upper bounds
        Scalar constr_viol_tol_ = 1e-4;          ///< Constraint violation tolerance
        Scalar bound_relax_factor_ = 1e-8;       ///< Factor for relaxing bounds
        IpTimingStatistics *timings_ = nullptr;  ///< Timing statistics

    public:
        // Setter methods for options
        void set_constr_viol_tol(const Scalar &value) { constr_viol_tol_ = value; }
        void set_bound_relax_factor(const Scalar &value) { bound_relax_factor_ = value; }

        void set_timing_statistics(IpTimingStatistics *timings) { timings_ = timings; }

        // Register options
        void register_options(OptionRegistry &registry);
    };

} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_nlp_orig_hpp__
