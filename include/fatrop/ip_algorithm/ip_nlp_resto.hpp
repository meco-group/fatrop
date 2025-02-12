//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_nlp_resto_hpp__
#define __fatrop_ip_algorithm_ip_nlp_resto_hpp__
#include "fatrop/common/fwd.hpp"
#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/nlp/dims.hpp"
#include "fatrop/nlp/fwd.hpp"
#include "fatrop/nlp/nlp.hpp"
#include <memory>
namespace fatrop
{
    /**
     * @brief Resto NLP formulation for interior point methods.
     *
     * This class wraps the original NLP problem and provides the interface
     * required for the interior point algorithm.
     *
     * @tparam ProblemType The type of optimization problem being solved.
     */
    template <typename ProblemType> class IpNlpResto : public Nlp<ProblemType>
    {
        typedef std::shared_ptr<Nlp<ProblemType>> NlpSp;

    public:
        /**
         * @brief Construct a new IpNlpResto object.
         *
         * @param nlp Shared pointer to the resto NLP problem.
         */
        IpNlpResto(const NlpSp &nlp);

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
        void set_xs_reference(const ProblemInfo<ProblemType> &info, const VecRealView &x_reference,
                              const VecRealView &s_reference);
        void set_zeta(const Scalar zeta) { zeta_ = zeta; };
        void set_rho(const Scalar rho) { rho_ = rho; };

    private:
        NlpSp nlp_; ///< Shared pointer to the original NLP problem
        NlpDims dims_;
        VecRealAllocated x_reference_;
        VecRealAllocated dr_;
        Scalar zeta_ = 1.;
        Scalar rho_ = 1.;

    public:
        // Register options
        void register_options(OptionRegistry &registry);
    };

} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_nlp_resto_hpp__
