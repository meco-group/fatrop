//
// Checks that the algorithm evaluates gradient, Jacobian and Hessian with one eval_derivatives
// call per accepted step, at the final multipliers, that an Nlp without an override gives the
// same solution, and that the restoration NLP's combined evaluation matches its individual ones.
//

#include "dense_test_problem.hpp"

#include "fatrop/common/options.hpp"
#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/nlp_dense.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/ip_alg_builder.hpp"
#include "fatrop/ip_algorithm/ip_algorithm.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_iterate.hpp"
#include "fatrop/ip_algorithm/ip_nlp_resto.hpp"
#include "fatrop/linear_algebra/matrix.hpp"
#include "fatrop/linear_algebra/vector.hpp"

#include <gtest/gtest.h>
#include <memory>
#include <vector>

using namespace fatrop;
using namespace fatrop::test;

namespace
{
    // Forwards to an inner Nlp and counts the derivative evaluations.
    class CountingNlp : public Nlp<DenseType>
    {
    public:
        explicit CountingNlp(std::shared_ptr<Nlp<DenseType>> inner) : inner_(std::move(inner)) {}

        const NlpDims &nlp_dims() const override { return inner_->nlp_dims(); }
        const ProblemDims<DenseType> &problem_dims() const override
        {
            return inner_->problem_dims();
        }
        Index eval_lag_hess(const ProblemInfo<DenseType> &info, const Scalar objective_scale,
                            const VecRealView &primal_x, const VecRealView &primal_s,
                            const VecRealView &lam, Hessian<DenseType> &hess) override
        {
            ++hess_calls;
            return inner_->eval_lag_hess(info, objective_scale, primal_x, primal_s, lam, hess);
        }
        Index eval_constr_jac(const ProblemInfo<DenseType> &info, const VecRealView &primal_x,
                              const VecRealView &primal_s, Jacobian<DenseType> &jac) override
        {
            ++jac_calls;
            return inner_->eval_constr_jac(info, primal_x, primal_s, jac);
        }
        Index eval_constraint_violation(const ProblemInfo<DenseType> &info,
                                        const VecRealView &primal_x, const VecRealView &primal_s,
                                        VecRealView &res) override
        {
            return inner_->eval_constraint_violation(info, primal_x, primal_s, res);
        }
        Index eval_objective_gradient(const ProblemInfo<DenseType> &info,
                                      const Scalar objective_scale, const VecRealView &primal_x,
                                      const VecRealView &primal_s, VecRealView &grad_x,
                                      VecRealView &grad_s) override
        {
            ++grad_calls;
            return inner_->eval_objective_gradient(info, objective_scale, primal_x, primal_s,
                                                   grad_x, grad_s);
        }
        Index eval_objective(const ProblemInfo<DenseType> &info, const Scalar objective_scale,
                             const VecRealView &primal_x, const VecRealView &primal_s,
                             Scalar &res) override
        {
            return inner_->eval_objective(info, objective_scale, primal_x, primal_s, res);
        }
        Index get_bounds(const ProblemInfo<DenseType> &info, VecRealView &lower_bounds,
                         VecRealView &upper_bounds) override
        {
            return inner_->get_bounds(info, lower_bounds, upper_bounds);
        }
        Index get_initial_primal(const ProblemInfo<DenseType> &info, VecRealView &primal_x) override
        {
            return inner_->get_initial_primal(info, primal_x);
        }
        void get_primal_damping(const ProblemInfo<DenseType> &info, VecRealView &damping) override
        {
            inner_->get_primal_damping(info, damping);
        }
        void apply_jacobian_s_transpose(const ProblemInfo<DenseType> &info,
                                        const VecRealView &multipliers, const Scalar alpha,
                                        const VecRealView &y, VecRealView &out) override
        {
            inner_->apply_jacobian_s_transpose(info, multipliers, alpha, y, out);
        }
        void callback(const IpData<DenseType> &ip_data) override { ++iterations; }

        Index hess_calls = 0;
        Index jac_calls = 0;
        Index grad_calls = 0;
        Index iterations = 0;

    protected:
        std::shared_ptr<Nlp<DenseType>> inner_;
    };

    // Additionally overrides eval_derivatives; forwards to the inner default implementation.
    class CombinedNlp : public CountingNlp
    {
    public:
        using CountingNlp::CountingNlp;
        Index eval_derivatives(const ProblemInfo<DenseType> &info, const Scalar objective_scale,
                               const VecRealView &primal_x, const VecRealView &primal_s,
                               const VecRealView &lam, VecRealView &grad_x, VecRealView &grad_s,
                               Jacobian<DenseType> &jac, Hessian<DenseType> &hess) override
        {
            ++combined_calls;
            last_lam.assign(lam.m(), 0.);
            for (Index i = 0; i < lam.m(); ++i)
                last_lam[i] = lam(i);
            return inner_->eval_derivatives(info, objective_scale, primal_x, primal_s, lam, grad_x,
                                            grad_s, jac, hess);
        }
        void callback(const IpData<DenseType> &ip_data) override
        {
            ++iterations;
            if (combined_calls == 0)
                return;
            // the derivatives used in this iteration were evaluated at the current multipliers
            const VecRealView &dual_eq =
                const_cast<IpData<DenseType> &>(ip_data).current_iterate().dual_eq();
            ASSERT_EQ(static_cast<Index>(last_lam.size()), dual_eq.m());
            for (Index i = 0; i < dual_eq.m(); ++i)
                EXPECT_EQ(last_lam[i], dual_eq(i));
        }

        Index combined_calls = 0;
        std::vector<Scalar> last_lam;
    };

    template <typename NlpType> std::shared_ptr<NlpType> solve(std::vector<Scalar> &solution)
    {
        auto nlp = std::make_shared<NlpType>(
            std::make_shared<NlpDense>(std::make_shared<DenseTestProblem>()));
        OptionRegistry options;
        IpAlgBuilder<DenseType> builder(nlp);
        auto alg = builder.with_options_registry(&options).build();
        EXPECT_EQ(static_cast<int>(alg->optimize()), static_cast<int>(IpSolverReturnFlag::Success));
        const VecRealView &x = alg->solution_primal();
        solution.assign(x.m(), 0.);
        for (Index i = 0; i < x.m(); ++i)
            solution[i] = x(i);
        return nlp;
    }
    void expect_mat_eq(const MatRealAllocated &ma, const MatRealAllocated &mb)
    {
        ASSERT_EQ(ma.m(), mb.m());
        ASSERT_EQ(ma.n(), mb.n());
        for (Index i = 0; i < ma.m(); ++i)
            for (Index j = 0; j < ma.n(); ++j)
                EXPECT_EQ(ma(i, j), mb(i, j));
    }
} // namespace

TEST(EvalDerivatives, OneCombinedCallPerAcceptedStep)
{
    std::vector<Scalar> x_default, x_combined;
    auto def = solve<CountingNlp>(x_default);
    auto comb = solve<CombinedNlp>(x_combined);

    // the default implementation reproduces the separate evaluations exactly
    ASSERT_GT(def->iterations, 0);
    EXPECT_EQ(def->iterations, comb->iterations);
    ASSERT_EQ(x_default.size(), x_combined.size());
    for (std::size_t i = 0; i < x_default.size(); ++i)
        EXPECT_EQ(x_default[i], x_combined[i]);

    // every accepted step (one per iteration) evaluates all derivatives in one call
    EXPECT_EQ(comb->combined_calls, comb->iterations);
    // only the initial iterate uses the individual functions
    EXPECT_EQ(comb->hess_calls, 1);
    EXPECT_EQ(comb->jac_calls, 1);
    EXPECT_EQ(comb->grad_calls, 1);
}

TEST(EvalDerivatives, RestoMatchesIndividualEvaluations)
{
    auto nlp = std::make_shared<IpNlpResto<DenseType>>(
        std::make_shared<NlpDense>(std::make_shared<DenseTestProblem>()));
    auto data = std::make_shared<IpData<DenseType>>(nlp);
    IpIterate<DenseType> &a = data->current_iterate();
    IpIterate<DenseType> &b = data->trial_iterate();
    const NlpDims &dims = nlp->nlp_dims();
    VecRealAllocated x(dims.number_of_variables), s(dims.number_of_ineq_constraints),
        lam(dims.number_of_eq_constraints);
    for (Index i = 0; i < x.m(); ++i)
        x(i) = 0.3 + 0.1 * i;
    for (Index i = 0; i < s.m(); ++i)
        s(i) = 0.5 + 0.2 * i;
    for (Index i = 0; i < lam.m(); ++i)
        lam(i) = 1.0 - 0.7 * i;
    for (IpIterate<DenseType> *it : {&a, &b})
    {
        it->set_primal_x(x);
        it->set_primal_s(s);
        it->set_dual_eq(lam);
    }

    a.eval_derivatives();
    for (Index i = 0; i < x.m(); ++i)
        EXPECT_EQ(a.obj_gradient_x()(i), b.obj_gradient_x()(i));
    for (Index i = 0; i < s.m(); ++i)
        EXPECT_EQ(a.obj_gradient_s()(i), b.obj_gradient_s()(i));
    expect_mat_eq(a.hessian().Hht, b.hessian().Hht);
    expect_mat_eq(a.jacobian().Gg_eqt, b.jacobian().Gg_eqt);
    expect_mat_eq(a.jacobian().Gg_ineqt, b.jacobian().Gg_ineqt);
}
