//
// Regression test: solve the same NLP via the DenseType pipeline and via the
// OcpType pipeline (with K=1, nu=0). The augmented-system solver for DenseType
// is a hand-derived simplification of the OcpType code path, so the iteration
// sequences must be bit-equivalent.
//

#include "dense_test_problem.hpp"

#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/nlp_dense.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/ip_algorithm/ip_alg_builder.hpp"
#include "fatrop/ip_algorithm/ip_algorithm.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_iterate.hpp"
#include "fatrop/common/options.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/nlp_ocp.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include "fatrop/ocp/type.hpp"

#include <gtest/gtest.h>
#include <memory>
#include <vector>

using namespace fatrop;
using namespace fatrop::test;

namespace
{

// One snapshot of the iterate state we care about for the regression check.
struct IterationSnapshot
{
    Scalar mu;
    std::vector<Scalar> primal_x;
    std::vector<Scalar> primal_s;
    std::vector<Scalar> dual_eq;
    std::vector<Scalar> dual_bounds_l;
    std::vector<Scalar> dual_bounds_u;
};

template <typename Vec> std::vector<Scalar> to_std(const Vec &v)
{
    std::vector<Scalar> out(v.m());
    for (Index i = 0; i < v.m(); ++i)
        out[i] = v(i);
    return out;
}

template <typename ProblemType>
IterationSnapshot snapshot(const IpData<ProblemType> &data)
{
    const auto &it = const_cast<IpData<ProblemType> &>(data).current_iterate();
    return {it.mu(),     to_std(it.primal_x()),      to_std(it.primal_s()),
            to_std(it.dual_eq()), to_std(it.dual_bounds_l()), to_std(it.dual_bounds_u())};
}

// Wraps an existing Nlp and captures a snapshot of every iterate the algorithm
// visits via the callback hook.
template <typename ProblemType> class LoggingNlp : public Nlp<ProblemType>
{
public:
    explicit LoggingNlp(std::shared_ptr<Nlp<ProblemType>> inner) : inner_(std::move(inner)) {}

    const NlpDims &nlp_dims() const override { return inner_->nlp_dims(); }
    const ProblemDims<ProblemType> &problem_dims() const override
    {
        return inner_->problem_dims();
    }
    Index eval_lag_hess(const ProblemInfo<ProblemType> &info, const Scalar objective_scale,
                        const VecRealView &primal_x, const VecRealView &primal_s,
                        const VecRealView &lam, Hessian<ProblemType> &hess) override
    {
        return inner_->eval_lag_hess(info, objective_scale, primal_x, primal_s, lam, hess);
    }
    Index eval_constr_jac(const ProblemInfo<ProblemType> &info, const VecRealView &primal_x,
                          const VecRealView &primal_s, Jacobian<ProblemType> &jac) override
    {
        return inner_->eval_constr_jac(info, primal_x, primal_s, jac);
    }
    Index eval_constraint_violation(const ProblemInfo<ProblemType> &info,
                                    const VecRealView &primal_x, const VecRealView &primal_s,
                                    VecRealView &res) override
    {
        return inner_->eval_constraint_violation(info, primal_x, primal_s, res);
    }
    Index eval_objective_gradient(const ProblemInfo<ProblemType> &info,
                                  const Scalar objective_scale, const VecRealView &primal_x,
                                  const VecRealView &primal_s, VecRealView &grad_x,
                                  VecRealView &grad_s) override
    {
        return inner_->eval_objective_gradient(info, objective_scale, primal_x, primal_s, grad_x,
                                               grad_s);
    }
    Index eval_objective(const ProblemInfo<ProblemType> &info, const Scalar objective_scale,
                         const VecRealView &primal_x, const VecRealView &primal_s,
                         Scalar &res) override
    {
        return inner_->eval_objective(info, objective_scale, primal_x, primal_s, res);
    }
    Index get_bounds(const ProblemInfo<ProblemType> &info, VecRealView &lower_bounds,
                     VecRealView &upper_bounds) override
    {
        return inner_->get_bounds(info, lower_bounds, upper_bounds);
    }
    Index get_initial_primal(const ProblemInfo<ProblemType> &info, VecRealView &primal_x) override
    {
        return inner_->get_initial_primal(info, primal_x);
    }
    void get_primal_damping(const ProblemInfo<ProblemType> &info, VecRealView &damping) override
    {
        inner_->get_primal_damping(info, damping);
    }
    void apply_jacobian_s_transpose(const ProblemInfo<ProblemType> &info,
                                    const VecRealView &multipliers, const Scalar alpha,
                                    const VecRealView &y, VecRealView &out) override
    {
        inner_->apply_jacobian_s_transpose(info, multipliers, alpha, y, out);
    }

    void callback(const IpData<ProblemType> &ip_data) override
    {
        log.push_back(snapshot(ip_data));
    }

    std::vector<IterationSnapshot> log;

private:
    std::shared_ptr<Nlp<ProblemType>> inner_;
};

} // namespace

TEST(DenseVsOcpRegression, IterationSequenceMatches)
{
    // ---- Dense pipeline ---------------------------------------------------
    auto dense_problem = std::make_shared<DenseTestProblem>();
    auto dense_inner = std::make_shared<NlpDense>(dense_problem);
    auto dense_nlp = std::make_shared<LoggingNlp<DenseType>>(dense_inner);
    OptionRegistry dense_options;
    IpAlgBuilder<DenseType> dense_builder(dense_nlp);
    auto dense_alg = dense_builder.with_options_registry(&dense_options).build();
    const IpSolverReturnFlag dense_ret = dense_alg->optimize();
    ASSERT_EQ(static_cast<int>(dense_ret), static_cast<int>(IpSolverReturnFlag::Success));

    // ---- OCP (K=1, nu=0) pipeline ----------------------------------------
    auto ocp_problem = std::make_shared<DenseAsOcpTestProblem>();
    auto ocp_inner = std::make_shared<NlpOcp>(ocp_problem);
    auto ocp_nlp = std::make_shared<LoggingNlp<OcpType>>(ocp_inner);
    OptionRegistry ocp_options;
    IpAlgBuilder<OcpType> ocp_builder(ocp_nlp);
    auto ocp_alg = ocp_builder.with_options_registry(&ocp_options).build();
    const IpSolverReturnFlag ocp_ret = ocp_alg->optimize();
    ASSERT_EQ(static_cast<int>(ocp_ret), static_cast<int>(IpSolverReturnFlag::Success));

    // ---- Compare per-iteration snapshots ---------------------------------
    ASSERT_EQ(dense_nlp->log.size(), ocp_nlp->log.size())
        << "iteration counts differ between the dense and OCP solvers";

    // The same arithmetic is performed in the same order, so equality is
    // bit-for-bit. We still use a small tolerance to be robust against
    // floating-point reorderings the compiler might introduce.
    constexpr Scalar kTol = 1e-14;
    for (std::size_t k = 0; k < dense_nlp->log.size(); ++k)
    {
        const IterationSnapshot &d = dense_nlp->log[k];
        const IterationSnapshot &o = ocp_nlp->log[k];
        EXPECT_NEAR(d.mu, o.mu, kTol) << "iter " << k << ": mu";
        ASSERT_EQ(d.primal_x.size(), o.primal_x.size());
        ASSERT_EQ(d.primal_s.size(), o.primal_s.size());
        ASSERT_EQ(d.dual_eq.size(), o.dual_eq.size());
        ASSERT_EQ(d.dual_bounds_l.size(), o.dual_bounds_l.size());
        ASSERT_EQ(d.dual_bounds_u.size(), o.dual_bounds_u.size());
        for (std::size_t i = 0; i < d.primal_x.size(); ++i)
            EXPECT_NEAR(d.primal_x[i], o.primal_x[i], kTol)
                << "iter " << k << ": primal_x[" << i << "]";
        for (std::size_t i = 0; i < d.primal_s.size(); ++i)
            EXPECT_NEAR(d.primal_s[i], o.primal_s[i], kTol)
                << "iter " << k << ": primal_s[" << i << "]";
        for (std::size_t i = 0; i < d.dual_eq.size(); ++i)
            EXPECT_NEAR(d.dual_eq[i], o.dual_eq[i], kTol)
                << "iter " << k << ": dual_eq[" << i << "]";
        for (std::size_t i = 0; i < d.dual_bounds_l.size(); ++i)
            EXPECT_NEAR(d.dual_bounds_l[i], o.dual_bounds_l[i], kTol)
                << "iter " << k << ": dual_bounds_l[" << i << "]";
        for (std::size_t i = 0; i < d.dual_bounds_u.size(); ++i)
            EXPECT_NEAR(d.dual_bounds_u[i], o.dual_bounds_u[i], kTol)
                << "iter " << k << ": dual_bounds_u[" << i << "]";
    }
}
