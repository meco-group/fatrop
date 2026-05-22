#ifndef __fatrop_unittest_dense_dense_test_problem_hpp__
#define __fatrop_unittest_dense_dense_test_problem_hpp__

#include "fatrop/dense/dense_abstract.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/ocp/ocp_abstract.hpp"
#include <limits>

// Shared dense NLP test problem.
//
//   variables : x in R^4
//   objective : 0.5 x' Q x + c' x, Q = diag(1, 2, 3, 4), c = (-1, -2, -3, -4)
//   eq        : x_0 + x_1 + x_2 + x_3 = 1
//               x_0 - x_3            = 0
//   ineq      : 0 <= x_0, x_1 <= 0.5
//
// Bounds on the (slack) inequalities: lower = 0, upper = +inf for x_0;
// lower = 0, upper = 0.5 for x_1. Two values to exercise both single-sided
// and double-sided bounds.

namespace fatrop
{
namespace test
{

// Convex quadratic objective + linear constraints. We model the math once and
// expose it through both the OcpAbstract (with K=1, nu=0) and DenseAbstract
// interfaces so the regression test can compare iteration sequences.
struct DenseTestProblemMath
{
    static constexpr Index nx = 4;
    static constexpr Index ng = 2;
    static constexpr Index ng_ineq = 2;

    // Q = diag(1, 2, 3, 4).
    static Scalar Q(Index i) { return Scalar(i + 1); }
    static Scalar c(Index i) { return -Scalar(i + 1); }

    static void eval_grad(const Scalar *objective_scale, const Scalar *x, Scalar *res)
    {
        const Scalar s = objective_scale[0];
        for (Index i = 0; i < nx; ++i)
            res[i] = s * (Q(i) * x[i] + c(i));
    }

    static void eval_f(const Scalar *objective_scale, const Scalar *x, Scalar *res)
    {
        Scalar v = 0.0;
        for (Index i = 0; i < nx; ++i)
            v += 0.5 * Q(i) * x[i] * x[i] + c(i) * x[i];
        *res = objective_scale[0] * v;
    }

    static void fill_hessian(const Scalar *objective_scale, MAT *Hh)
    {
        const Scalar s = objective_scale[0];
        blasfeo_gese_wrap(Hh->m, Hh->n, 0.0, Hh, 0, 0);
        for (Index i = 0; i < nx; ++i)
            blasfeo_matel_wrap(Hh, i, i) = s * Q(i);
    }

    // Equality Jacobian (nx + 1) x ng. Top nx rows: ∇g_eq^T. Last row: g_eq(x).
    static void fill_Ggt(const Scalar *x, MAT *Ggt)
    {
        blasfeo_gese_wrap(Ggt->m, Ggt->n, 0.0, Ggt, 0, 0);
        // Row 0: x_0 + x_1 + x_2 + x_3 - 1
        for (Index i = 0; i < nx; ++i)
            blasfeo_matel_wrap(Ggt, i, 0) = 1.0;
        blasfeo_matel_wrap(Ggt, nx, 0) = x[0] + x[1] + x[2] + x[3] - 1.0;
        // Row 1: x_0 - x_3
        blasfeo_matel_wrap(Ggt, 0, 1) = 1.0;
        blasfeo_matel_wrap(Ggt, 3, 1) = -1.0;
        blasfeo_matel_wrap(Ggt, nx, 1) = x[0] - x[3];
    }

    static void eval_g(const Scalar *x, Scalar *res)
    {
        res[0] = x[0] + x[1] + x[2] + x[3] - 1.0;
        res[1] = x[0] - x[3];
    }

    // Inequality Jacobian (nx + 1) x ng_ineq.
    static void fill_Ggt_ineq(const Scalar *x, MAT *Ggt)
    {
        blasfeo_gese_wrap(Ggt->m, Ggt->n, 0.0, Ggt, 0, 0);
        blasfeo_matel_wrap(Ggt, 0, 0) = 1.0; // g_ineq_0 = x_0
        blasfeo_matel_wrap(Ggt, 1, 1) = 1.0; // g_ineq_1 = x_1
        blasfeo_matel_wrap(Ggt, nx, 0) = x[0];
        blasfeo_matel_wrap(Ggt, nx, 1) = x[1];
    }

    static void eval_gineq(const Scalar *x, Scalar *res)
    {
        res[0] = x[0];
        res[1] = x[1];
    }

    static void get_bounds(Scalar *lower, Scalar *upper)
    {
        lower[0] = 0.0;
        upper[0] = std::numeric_limits<Scalar>::infinity(); // single-sided
        lower[1] = 0.0;
        upper[1] = 0.5; // double-sided
    }

    static void get_initial(Scalar *x)
    {
        x[0] = 0.1;
        x[1] = 0.2;
        x[2] = 0.3;
        x[3] = 0.4;
    }
};

class DenseTestProblem : public DenseAbstract
{
public:
    Index get_nx() const override { return DenseTestProblemMath::nx; }
    Index get_ng() const override { return DenseTestProblemMath::ng; }
    Index get_ng_ineq() const override { return DenseTestProblemMath::ng_ineq; }

    Index eval_Hh(const Scalar *objective_scale, const Scalar *x, const Scalar * /*lam*/,
                  MAT *res) override
    {
        DenseTestProblemMath::fill_hessian(objective_scale, res);
        Scalar grad[DenseTestProblemMath::nx];
        DenseTestProblemMath::eval_grad(objective_scale, x, grad);
        for (Index i = 0; i < DenseTestProblemMath::nx; ++i)
            blasfeo_matel_wrap(res, DenseTestProblemMath::nx, i) = grad[i];
        return 0;
    }

    Index eval_Ggt(const Scalar *x, MAT *res) override
    {
        DenseTestProblemMath::fill_Ggt(x, res);
        return 0;
    }

    Index eval_Ggt_ineq(const Scalar *x, MAT *res) override
    {
        DenseTestProblemMath::fill_Ggt_ineq(x, res);
        return 0;
    }

    Index eval_g(const Scalar *x, Scalar *res) override
    {
        DenseTestProblemMath::eval_g(x, res);
        return 0;
    }
    Index eval_gineq(const Scalar *x, Scalar *res) override
    {
        DenseTestProblemMath::eval_gineq(x, res);
        return 0;
    }
    Index eval_grad(const Scalar *objective_scale, const Scalar *x, Scalar *res) override
    {
        DenseTestProblemMath::eval_grad(objective_scale, x, res);
        return 0;
    }
    Index eval_f(const Scalar *objective_scale, const Scalar *x, Scalar *res) override
    {
        DenseTestProblemMath::eval_f(objective_scale, x, res);
        return 0;
    }
    Index get_bounds(Scalar *lower, Scalar *upper) const override
    {
        DenseTestProblemMath::get_bounds(lower, upper);
        return 0;
    }
    Index get_initial(Scalar *x) const override
    {
        DenseTestProblemMath::get_initial(x);
        return 0;
    }
};

// Same math, exposed as a K=1, nu=0 OCP. Used by the regression test to verify
// that the DenseType solver produces bit-equivalent iterations to the OCP solver.
class DenseAsOcpTestProblem : public OcpAbstract
{
public:
    Index get_nx(const Index) const override { return DenseTestProblemMath::nx; }
    Index get_nu(const Index) const override { return 0; }
    Index get_ng(const Index) const override { return DenseTestProblemMath::ng; }
    Index get_ng_ineq(const Index) const override { return DenseTestProblemMath::ng_ineq; }
    Index get_horizon_length() const override { return 1; }

    // K = 1 so there is no dynamics block.
    Index eval_BAbt(const Scalar *, const Scalar *, const Scalar *, MAT *, const Index) override
    {
        return 0;
    }
    Index eval_b(const Scalar *, const Scalar *, const Scalar *, Scalar *, const Index) override
    {
        return 0;
    }

    Index eval_RSQrqt(const Scalar *objective_scale, const Scalar * /*inputs_k*/,
                      const Scalar *states_k, const Scalar * /*lam_dyn_k*/,
                      const Scalar * /*lam_eq_k*/, const Scalar * /*lam_eq_ineq_k*/, MAT *res,
                      const Index) override
    {
        DenseTestProblemMath::fill_hessian(objective_scale, res);
        Scalar grad[DenseTestProblemMath::nx];
        DenseTestProblemMath::eval_grad(objective_scale, states_k, grad);
        for (Index i = 0; i < DenseTestProblemMath::nx; ++i)
            blasfeo_matel_wrap(res, DenseTestProblemMath::nx, i) = grad[i];
        return 0;
    }

    Index eval_Ggt(const Scalar * /*inputs_k*/, const Scalar *states_k, MAT *res,
                   const Index) override
    {
        DenseTestProblemMath::fill_Ggt(states_k, res);
        return 0;
    }

    Index eval_Ggt_ineq(const Scalar * /*inputs_k*/, const Scalar *states_k, MAT *res,
                        const Index) override
    {
        DenseTestProblemMath::fill_Ggt_ineq(states_k, res);
        return 0;
    }

    Index eval_g(const Scalar * /*inputs_k*/, const Scalar *states_k, Scalar *res,
                 const Index) override
    {
        DenseTestProblemMath::eval_g(states_k, res);
        return 0;
    }

    Index eval_gineq(const Scalar * /*inputs_k*/, const Scalar *states_k, Scalar *res,
                     const Index) override
    {
        DenseTestProblemMath::eval_gineq(states_k, res);
        return 0;
    }

    Index eval_rq(const Scalar *objective_scale, const Scalar * /*inputs_k*/,
                  const Scalar *states_k, Scalar *res, const Index) override
    {
        DenseTestProblemMath::eval_grad(objective_scale, states_k, res);
        return 0;
    }

    Index eval_L(const Scalar *objective_scale, const Scalar * /*inputs_k*/,
                 const Scalar *states_k, Scalar *res, const Index) override
    {
        DenseTestProblemMath::eval_f(objective_scale, states_k, res);
        return 0;
    }

    Index get_bounds(Scalar *lower, Scalar *upper, const Index) const override
    {
        DenseTestProblemMath::get_bounds(lower, upper);
        return 0;
    }

    Index get_initial_xk(Scalar *xk, const Index) const override
    {
        DenseTestProblemMath::get_initial(xk);
        return 0;
    }
    Index get_initial_uk(Scalar *, const Index) const override { return 0; }
};

} // namespace test
} // namespace fatrop

#endif // __fatrop_unittest_dense_dense_test_problem_hpp__
