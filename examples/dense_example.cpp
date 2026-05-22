//
// Dense NLP example: a small constrained quadratic problem solved via the
// new DenseAbstract / NlpDense / DenseType pipeline.
//
//   min   f(x) = 0.5 * (x1^2 + 2*x2^2 + 3*x3^2) - x1 - 2*x2 - 3*x3
//   s.t.  x1 + x2 + x3 = 1                  (equality)
//         0 <= x1, x2, x3 <= 1              (3 box-style inequalities)
//
// The objective is strictly convex, so the optimum is unique. The unconstrained
// stationarity equations alone would give x1 < 0, so the lower bound on x1 is
// active at the optimum. With x1 = 0 the eq constraint and the remaining
// stationarity conditions give x* = (0, 2/5, 3/5), eq multiplier λ* = 6/5,
// and bound multiplier μ_1 = 1/5 (verified by hand).
//

#include <cmath>
#include <iostream>
#include <memory>
#include <fatrop/fatrop.hpp>

using namespace fatrop;

class DenseQpExample : public DenseAbstract
{
public:
    Index get_nx() const override { return 3; }
    Index get_ng() const override { return 1; }
    Index get_ng_ineq() const override { return 3; }

    // Hht is (nx + 1) x nx. Top nx-by-nx is the Hessian, last row is the gradient.
    Index eval_Hh(const Scalar *objective_scale, const Scalar *x, const Scalar * /*lam*/,
                  MAT *res) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        const Scalar s = objective_scale[0];
        // Hessian H = diag(1, 2, 3) (objective is quadratic, lam doesn't enter for linear constraints).
        blasfeo_matel_wrap(res, 0, 0) = s * 1.0;
        blasfeo_matel_wrap(res, 1, 1) = s * 2.0;
        blasfeo_matel_wrap(res, 2, 2) = s * 3.0;
        // Gradient row: ∇f(x) = (x1 - 1, 2*x2 - 2, 3*x3 - 3).
        blasfeo_matel_wrap(res, 3, 0) = s * (x[0] - 1.0);
        blasfeo_matel_wrap(res, 3, 1) = s * (2.0 * x[1] - 2.0);
        blasfeo_matel_wrap(res, 3, 2) = s * (3.0 * x[2] - 3.0);
        return 0;
    }

    // Equality Jacobian: ∇(x1 + x2 + x3) = [1, 1, 1].
    Index eval_Ggt(const Scalar *x, MAT *res) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        blasfeo_matel_wrap(res, 0, 0) = 1.0;
        blasfeo_matel_wrap(res, 1, 0) = 1.0;
        blasfeo_matel_wrap(res, 2, 0) = 1.0;
        // last row holds the equality residual g(x) = x1 + x2 + x3 - 1
        blasfeo_matel_wrap(res, 3, 0) = x[0] + x[1] + x[2] - 1.0;
        return 0;
    }

    // Inequality Jacobian: ∇(x_i) = e_i, so the Jacobian is identity.
    Index eval_Ggt_ineq(const Scalar *x, MAT *res) override
    {
        blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
        blasfeo_matel_wrap(res, 0, 0) = 1.0;
        blasfeo_matel_wrap(res, 1, 1) = 1.0;
        blasfeo_matel_wrap(res, 2, 2) = 1.0;
        // last row: g_ineq(x) = (x1, x2, x3)
        blasfeo_matel_wrap(res, 3, 0) = x[0];
        blasfeo_matel_wrap(res, 3, 1) = x[1];
        blasfeo_matel_wrap(res, 3, 2) = x[2];
        return 0;
    }

    Index eval_g(const Scalar *x, Scalar *res) override
    {
        res[0] = x[0] + x[1] + x[2] - 1.0;
        return 0;
    }

    Index eval_gineq(const Scalar *x, Scalar *res) override
    {
        res[0] = x[0];
        res[1] = x[1];
        res[2] = x[2];
        return 0;
    }

    Index eval_grad(const Scalar *objective_scale, const Scalar *x, Scalar *res) override
    {
        const Scalar s = objective_scale[0];
        res[0] = s * (x[0] - 1.0);
        res[1] = s * (2.0 * x[1] - 2.0);
        res[2] = s * (3.0 * x[2] - 3.0);
        return 0;
    }

    Index eval_f(const Scalar *objective_scale, const Scalar *x, Scalar *res) override
    {
        *res = objective_scale[0] * (0.5 * (x[0] * x[0] + 2.0 * x[1] * x[1] + 3.0 * x[2] * x[2]) -
                                     x[0] - 2.0 * x[1] - 3.0 * x[2]);
        return 0;
    }

    Index get_bounds(Scalar *lower, Scalar *upper) const override
    {
        lower[0] = 0.0;
        upper[0] = 1.0;
        lower[1] = 0.0;
        upper[1] = 1.0;
        lower[2] = 0.0;
        upper[2] = 1.0;
        return 0;
    }

    Index get_initial(Scalar *x) const override
    {
        x[0] = 0.5;
        x[1] = 0.5;
        x[2] = 0.5;
        return 0;
    }
};

namespace
{
// Pretty-print and report whether x matches the closed-form optimum.
int report(const char *label, IpSolverReturnFlag ret, const VecRealView &x)
{
    const Scalar x_star[3] = {0.0, 2.0 / 5.0, 3.0 / 5.0}; // closed form
    Scalar err = 0.0;
    for (int i = 0; i < 3; ++i)
        err = std::max(err, std::abs(x(i) - x_star[i]));
    std::cout << "[" << label << "] return flag = " << int(ret)
              << "  x* = (" << x(0) << ", " << x(1) << ", " << x(2) << ")"
              << "  max abs error vs closed form = " << err << "\n";
    return (ret == IpSolverReturnFlag::Success && err < 1e-6) ? 0 : 1;
}
} // namespace

int main()
{
    auto problem = std::make_shared<DenseQpExample>();

    // ----------------------------- NLP (fatrop) -----------------------------
    OptionRegistry nlp_options;
    IpAlgBuilder<DenseType> nlp_builder(std::make_shared<NlpDense>(problem));
    auto nlp_alg = nlp_builder.with_options_registry(&nlp_options).build();

    Timer nlp_timer;
    nlp_timer.start();
    IpSolverReturnFlag nlp_ret = nlp_alg->optimize();
    const auto nlp_elapsed = nlp_timer.stop();
    auto nlp_data = nlp_builder.get_ipdata();
    std::cout << "\n=== Fatrop NLP (interior point) ===\n";
    std::cout << "Elapsed time: " << nlp_elapsed << "\n";
    std::cout << nlp_data->timing_statistics() << "\n";
    int nlp_status = report("NLP", nlp_ret, nlp_data->current_iterate().primal_x());

    // ----------------------------- Mehrotra QP ------------------------------
    // The objective is quadratic and the constraints are linear, so the same
    // problem is also a strictly convex QP. The Mehrotra predictor-corrector
    // solver should converge in a handful of iterations.
    OptionRegistry qp_options;
    MehrotraQpBuilder<DenseType> qp_builder(std::make_shared<NlpDense>(problem));
    auto qp_alg = qp_builder.with_options_registry(&qp_options).build();

    Timer qp_timer;
    qp_timer.start();
    IpSolverReturnFlag qp_ret = qp_alg->optimize();
    const auto qp_elapsed = qp_timer.stop();
    auto qp_data = qp_builder.get_ipdata();
    std::cout << "\n=== Mehrotra QP (predictor-corrector) ===\n";
    std::cout << "Elapsed time: " << qp_elapsed << "\n";
    std::cout << qp_data->timing_statistics() << "\n";
    int qp_status = report("QP ", qp_ret, qp_data->current_iterate().primal_x());

    return (nlp_status == 0 && qp_status == 0) ? 0 : 1;
}
