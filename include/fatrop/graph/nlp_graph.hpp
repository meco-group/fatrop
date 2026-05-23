//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_graph_nlp_graph_hpp__
#define __fatrop_graph_nlp_graph_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/graph/dims.hpp"
#include "fatrop/graph/fwd.hpp"
#include "fatrop/graph/graph_abstract.hpp"
#include "fatrop/graph/problem_info.hpp"
#include "fatrop/graph/problem_type.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/nlp/dims.hpp"
#include "fatrop/nlp/fwd.hpp"
#include "fatrop/nlp/nlp.hpp"

#include <memory>

namespace fatrop
{
    /**
     * @brief Adapter from a @c GraphAbstract user problem to fatrop's
     *        @c Nlp<GraphProblem>.
     *
     * Analogue of @c NlpOcp / @c NlpDense; it walks the per-block evaluation
     * hooks of @c GraphAbstract and writes the results into the block-sparse
     * Hessian / gradient and per-block transposed inequality Jacobian.
     */
    class NlpGraph : public Nlp<GraphProblem>
    {
        typedef std::shared_ptr<GraphAbstract> GraphAbstractSp;
        typedef ProblemInfo<GraphProblem> GraphInfo;

    public:
        NlpGraph(const GraphAbstractSp &graph);

        const NlpDims &nlp_dims() const override { return nlp_dims_; }
        const ProblemDims<GraphProblem> &problem_dims() const override { return dims_; }

        Index eval_lag_hess(const GraphInfo &info, const Scalar objective_scale,
                            const VecRealView &primal_x, const VecRealView &primal_s,
                            const VecRealView &lam, Hessian<GraphProblem> &hess) override;
        Index eval_constr_jac(const GraphInfo &info, const VecRealView &primal_x,
                              const VecRealView &primal_s, Jacobian<GraphProblem> &jac) override;
        Index eval_constraint_violation(const GraphInfo &info, const VecRealView &primal_x,
                                        const VecRealView &primal_s, VecRealView &res) override;
        Index eval_objective_gradient(const GraphInfo &info, const Scalar objective_scale,
                                      const VecRealView &primal_x, const VecRealView &primal_s,
                                      VecRealView &grad_x, VecRealView &grad_s) override;
        Index eval_objective(const GraphInfo &info, const Scalar objective_scale,
                             const VecRealView &primal_x, const VecRealView &primal_s,
                             Scalar &res) override;
        Index get_bounds(const GraphInfo &info, VecRealView &lower_bounds,
                         VecRealView &upper_bounds) override;
        Index get_initial_primal(const GraphInfo &info, VecRealView &primal_x) override;
        void get_primal_damping(const GraphInfo &info, VecRealView &damping) override;
        void apply_jacobian_s_transpose(const GraphInfo &info, const VecRealView &multipliers,
                                        const Scalar alpha, const VecRealView &y,
                                        VecRealView &out) override;
        void apply_retraction(const GraphInfo &info, const VecRealView &primal_x,
                              const VecRealView &delta_primal_x, const Scalar alpha,
                              VecRealView &primal_x_next) override;

        /**
         * @brief Set a constant Levenberg-Marquardt-style primal damping.
         *
         * Adds @c lambda * I to the Hessian on every IP iteration. This is
         * useful for graph problems where the Gauss-Newton Hessian
         * @c J^T J has near-zero eigenvalues in poorly-observed directions
         * (e.g. landmark depth from nearby cameras): without damping the
         * Newton step overshoots in those directions and the line search
         * spends many backtracks recovering. A small @c lambda (typically
         * 1e-3 to 1e0 times the smallest meaningful Hessian eigenvalue)
         * regularises the step without biasing the optimum significantly.
         *
         * Default is 0 (no damping).
         */
        void set_primal_damping_lambda(Scalar lambda) { lm_lambda_ = lambda; }

    private:
        GraphAbstractSp graph_;
        ProblemDims<GraphProblem> dims_;
        NlpDims nlp_dims_;
        Scalar lm_lambda_ = 0.0;
    };
} // namespace fatrop

#endif // __fatrop_graph_nlp_graph_hpp__
