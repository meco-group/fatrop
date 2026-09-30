#ifndef FATROP_BENCHMARK_OCP_DERIVATIVE_CHECKER_HPP
#define FATROP_BENCHMARK_OCP_DERIVATIVE_CHECKER_HPP

#include <fatrop/linear_algebra/matrix.hpp>
#include <fatrop/ocp/ocp_abstract.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace fatrop_benchmark
{
struct DerivativeReport
{
    int checks = 0, failures = 0;
    double max_scaled_error = 0.;
    std::string worst;

    void compare(const std::string &place, double analytic, double finite_difference)
    {
        ++checks;
        const double error = std::abs(analytic - finite_difference);
        const double scale = 1. + std::max(std::abs(analytic), std::abs(finite_difference));
        const double scaled = error / scale;
        if (scaled > max_scaled_error)
        {
            max_scaled_error = scaled;
            worst = place + " analytic=" + std::to_string(analytic)
                  + " fd=" + std::to_string(finite_difference);
        }
        if (!std::isfinite(scaled) || error > 2e-5 * scale) ++failures;
    }
};

// Compares OCP derivatives with finite differences at each stage.
class OcpDerivativeChecker
{
public:
    explicit OcpDerivativeChecker(fatrop::OcpAbstract &problem) : problem_(problem) {}

    DerivativeReport check_all_initial_points()
    {
        DerivativeReport report;
        for (fatrop::Index k = 0; k < problem_.get_horizon_length(); ++k)
        {
            const fatrop::Index nu = problem_.get_nu(k), nx = problem_.get_nx(k);
            if (nu != problem_.get_nu_tangent(k) || nx != problem_.get_nx_tangent(k))
                throw std::invalid_argument("finite differences require Euclidean OCP coordinates");
            std::vector<fatrop::Scalar> u(nu), x(nx), next(
                k + 1 < problem_.get_horizon_length() ? problem_.get_nx(k + 1) : 0);
            require(problem_.get_initial_uk(u.data(), k));
            require(problem_.get_initial_xk(x.data(), k));
            if (!next.empty()) require(problem_.get_initial_xk(next.data(), k + 1));
            // Move away from zeros, where odd nonlinearities can hide errors.
            for (fatrop::Index i = 0; i < nu; ++i)
                u[i] += 0.019 * std::cos(0.7 * (i + k + 1));
            for (fatrop::Index i = 0; i < nx; ++i)
                x[i] += 0.031 * std::sin(0.8 * (i + k + 1));
            check_stage(k, u, x, next, report);
        }
        return report;
    }

    void check_stage(fatrop::Index k, const std::vector<fatrop::Scalar> &u,
                     const std::vector<fatrop::Scalar> &x,
                     const std::vector<fatrop::Scalar> &next,
                     DerivativeReport &report)
    {
        const fatrop::Index nu = problem_.get_nu(k), nx = problem_.get_nx(k);
        const fatrop::Index nz = nu + nx;
        const fatrop::Index nd = k + 1 < problem_.get_horizon_length() ? problem_.get_nx(k + 1) : 0;
        const fatrop::Index ng = problem_.get_ng(k), ni = problem_.get_ng_ineq(k);
        if (fatrop::Index(u.size()) != nu || fatrop::Index(x.size()) != nx || fatrop::Index(next.size()) != nd)
            throw std::invalid_argument("finite-difference stage vectors have the wrong dimensions");
        if (nu != problem_.get_nu_tangent(k) || nx != problem_.get_nx_tangent(k))
            throw std::invalid_argument("finite differences require Euclidean OCP coordinates");

        std::vector<fatrop::Scalar> z(nz), plus(nz), minus(nz);
        std::copy(u.begin(), u.end(), z.begin());
        std::copy(x.begin(), x.end(), z.begin() + nu);
        auto ux = [&](const std::vector<fatrop::Scalar> &v) {
            return std::pair<const fatrop::Scalar *, const fatrop::Scalar *>(v.data(), v.data() + nu);
        };
        auto check_matrix = [&](const char *name, fatrop::Index count,
                                auto jacobian, auto residual) {
            if (count == 0) return;
            fatrop::MatRealAllocated analytic(nz + 1, count);
            analytic = 0.;
            auto [up, xp] = ux(z);
            require(jacobian(up, xp, &analytic.mat()));
            std::vector<fatrop::Scalar> high(count), low(count);
            for (fatrop::Index j = 0; j < nz; ++j)
            {
                plus = minus = z;
                const double step = 1e-6 * std::max(1., std::abs(z[j]));
                plus[j] += step; minus[j] -= step;
                auto [u_high, x_high] = ux(plus);
                auto [u_low, x_low] = ux(minus);
                require(residual(u_high, x_high, high.data()));
                require(residual(u_low, x_low, low.data()));
                for (fatrop::Index i = 0; i < count; ++i)
                    report.compare(label(name, k, j, i), analytic(j, i),
                                   (high[i] - low[i]) / (2. * step));
            }
        };
        check_matrix("dynamics Jacobian", nd,
            [&](const auto *up, const auto *xp, auto *out) {
                return problem_.eval_BAbt(next.data(), up, xp, out, k);
            },
            [&](const auto *up, const auto *xp, auto *out) {
                return problem_.eval_b(next.data(), up, xp, out, k);
            });
        check_matrix("equality Jacobian", ng,
            [&](const auto *up, const auto *xp, auto *out) {
                return problem_.eval_Ggt(up, xp, out, k);
            },
            [&](const auto *up, const auto *xp, auto *out) {
                return problem_.eval_g(up, xp, out, k);
            });
        check_matrix("inequality Jacobian", ni,
            [&](const auto *up, const auto *xp, auto *out) {
                return problem_.eval_Ggt_ineq(up, xp, out, k);
            },
            [&](const auto *up, const auto *xp, auto *out) {
                return problem_.eval_gineq(up, xp, out, k);
            });

        const fatrop::Scalar scale = 1.3;
        std::vector<fatrop::Scalar> gradient(nz), high_gradient(nz), low_gradient(nz);
        auto [up, xp] = ux(z);
        require(problem_.eval_rq(&scale, up, xp, gradient.data(), k));
        for (fatrop::Index j = 0; j < nz; ++j)
        {
            plus = minus = z;
            const double step = 1e-6 * std::max(1., std::abs(z[j]));
            plus[j] += step; minus[j] -= step;
            fatrop::Scalar high = 0., low = 0.;
            auto [u_high, x_high] = ux(plus);
            auto [u_low, x_low] = ux(minus);
            require(problem_.eval_L(&scale, u_high, x_high, &high, k));
            require(problem_.eval_L(&scale, u_low, x_low, &low, k));
            report.compare(label("objective gradient", k, j, 0), gradient[j],
                           (high - low) / (2. * step));
        }

        std::vector<fatrop::Scalar> lambda_dyn(nd), lambda_eq(ng), lambda_ineq(ni);
        for (fatrop::Index i = 0; i < nd; ++i) lambda_dyn[i] = 0.17 * (i + 1);
        for (fatrop::Index i = 0; i < ng; ++i) lambda_eq[i] = -0.11 * (i + 1);
        for (fatrop::Index i = 0; i < ni; ++i) lambda_ineq[i] = 0.09 * (i + 1);
        fatrop::MatRealAllocated hessian(nz + 1, nz);
        hessian = 0.;
        require(problem_.eval_RSQrqt(&scale, up, xp, lambda_dyn.data(),
                                     lambda_eq.data(), lambda_ineq.data(), &hessian.mat(), k));
        auto lag_gradient = [&](const std::vector<fatrop::Scalar> &point,
                                std::vector<fatrop::Scalar> &out) {
            auto [input, state] = ux(point);
            require(problem_.eval_rq(&scale, input, state, out.data(), k));
            auto add = [&](fatrop::Index count, const std::vector<fatrop::Scalar> &lambda,
                           auto evaluate) {
                if (!count) return;
                fatrop::MatRealAllocated jac(nz + 1, count);
                jac = 0.;
                require(evaluate(input, state, &jac.mat()));
                for (fatrop::Index j = 0; j < nz; ++j)
                    for (fatrop::Index i = 0; i < count; ++i)
                        out[j] += jac(j, i) * lambda[i];
            };
            add(nd, lambda_dyn, [&](const auto *a, const auto *b, auto *mat) {
                return problem_.eval_BAbt(next.data(), a, b, mat, k);
            });
            add(ng, lambda_eq, [&](const auto *a, const auto *b, auto *mat) {
                return problem_.eval_Ggt(a, b, mat, k);
            });
            add(ni, lambda_ineq, [&](const auto *a, const auto *b, auto *mat) {
                return problem_.eval_Ggt_ineq(a, b, mat, k);
            });
        };
        for (fatrop::Index j = 0; j < nz; ++j)
        {
            plus = minus = z;
            const double step = 1e-6 * std::max(1., std::abs(z[j]));
            plus[j] += step; minus[j] -= step;
            lag_gradient(plus, high_gradient);
            lag_gradient(minus, low_gradient);
            for (fatrop::Index i = j; i < nz; ++i)
                report.compare(label("Lagrangian Hessian", k, i, j), hessian(i, j),
                               (high_gradient[i] - low_gradient[i]) / (2. * step));
        }
    }

private:
    fatrop::OcpAbstract &problem_;
    static void require(fatrop::Index status)
    { if (status) throw std::runtime_error("OCP derivative callback failed"); }
    static std::string label(const char *kind, fatrop::Index stage,
                             fatrop::Index row, fatrop::Index column)
    {
        return std::string(kind) + " stage=" + std::to_string(stage)
             + " row=" + std::to_string(row) + " col=" + std::to_string(column);
    }
};
}
#endif
