//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_dense_nlp_dense_hpp__
#define __fatrop_dense_nlp_dense_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/dense/dense_abstract.hpp"
#include "fatrop/dense/dims.hpp"
#include "fatrop/dense/fwd.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/nlp/dims.hpp"
#include "fatrop/nlp/fwd.hpp"
#include "fatrop/nlp/nlp.hpp"

#include <memory>

namespace fatrop
{
    /**
     * @brief Adapter from a @c DenseAbstract user problem to fatrop's @c Nlp<DenseType>.
     *
     * This is the dense analogue of @c NlpOcp; same role, much simpler implementation
     * because there is no time/stage structure to unfold.
     */
    class NlpDense : public Nlp<DenseType>
    {
        typedef std::shared_ptr<DenseAbstract> DenseAbstractSp;
        typedef ProblemInfo<DenseType> DenseInfo;

    public:
        NlpDense(const DenseAbstractSp &dense);

        const NlpDims &nlp_dims() const override { return nlp_dims_; }
        const ProblemDims<DenseType> &problem_dims() const override { return dense_dims_; }

        Index eval_lag_hess(const DenseInfo &info, const Scalar objective_scale,
                            const VecRealView &primal_x, const VecRealView &primal_s,
                            const VecRealView &lam, Hessian<DenseType> &hess) override;
        Index eval_constr_jac(const DenseInfo &info, const VecRealView &primal_x,
                              const VecRealView &primal_s, Jacobian<DenseType> &jac) override;
        Index eval_constraint_violation(const DenseInfo &info, const VecRealView &primal_x,
                                        const VecRealView &primal_s, VecRealView &res) override;
        Index eval_objective_gradient(const DenseInfo &info, const Scalar objective_scale,
                                      const VecRealView &primal_x, const VecRealView &primal_s,
                                      VecRealView &grad_x, VecRealView &grad_s) override;
        Index eval_objective(const DenseInfo &info, const Scalar objective_scale,
                             const VecRealView &primal_x, const VecRealView &primal_s,
                             Scalar &res) override;
        Index get_bounds(const DenseInfo &info, VecRealView &lower_bounds,
                         VecRealView &upper_bounds) override;
        Index get_initial_primal(const DenseInfo &info, VecRealView &primal_x) override;
        void get_primal_damping(const DenseInfo &info, VecRealView &damping) override;
        void apply_jacobian_s_transpose(const DenseInfo &info, const VecRealView &multipliers,
                                        const Scalar alpha, const VecRealView &y,
                                        VecRealView &out) override;
        void apply_retraction(const DenseInfo &info, const VecRealView &primal_x,
                              const VecRealView &delta_primal_x, const Scalar alpha,
                              VecRealView &primal_x_next) override;
        void apply_dual_eq_transformation(const DenseInfo &info, const VecRealView &primal_x,
                                          const VecRealView &dual_eq_in,
                                          VecRealView &dual_eq_out) override;

    private:
        DenseAbstractSp dense_;
        ProblemDims<DenseType> dense_dims_;
        NlpDims nlp_dims_;
    };
} // namespace fatrop

#endif // __fatrop_dense_nlp_dense_hpp__
