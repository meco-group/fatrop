#include "../ocp_derivative_checker.hpp"
#include "../problems/problem_catalog.hpp"
#include "compiled_derivatives.hpp"
#include <iostream>
#include <cmath>
using namespace fatrop;
using namespace fatrop_benchmark;

template<class P> bool check()
{
    P reference;
    CompiledProblem<P> compiled;
    OcpDerivativeChecker checker(compiled);
    DerivativeReport fd;
    double worst = 0.0;
    for (int k = 0; k < P::stages; ++k)
    {
        const int nu = reference.get_nu(k), nx = P::states, n = nu + nx;
        const int ni = reference.get_ng_ineq(k), ng = reference.get_ng(k);
        const int nd = k + 1 < P::stages ? nx : 0;
        for (int sample = 0; sample < 3; ++sample)
        {
            std::vector<double> u(nu), x(nx), next(nd), ld(nd), li(ni), lg(ng);
            reference.get_initial_uk(u.data(), k);
            reference.get_initial_xk(x.data(), k);
            if (nd) reference.get_initial_xk(next.data(), k + 1);
            for (int j = 0; j < nu; ++j) u[j] += 0.13 * std::sin(0.8 * (j + k + sample + 1));
            for (int j = 0; j < nx; ++j) x[j] += 0.17 * std::cos(0.6 * (j + k + sample + 1));
            for (int j = 0; j < nd; ++j) ld[j] = std::sin(j + k + sample + 1.0);
            for (int j = 0; j < ni; ++j) li[j] = std::cos(j + k + sample + 1.0);
            const double scale = sample == 0 ? 0.0 : sample == 1 ? 1.3 : -0.7;
            MatRealAllocated a(n + 1, n), b(n + 1, n);
            a = 0.; b = 0.;
            reference.eval_RSQrqt(&scale, u.data(), x.data(), ld.data(), lg.data(), li.data(), &a.mat(), k);
            compiled.eval_RSQrqt(&scale, u.data(), x.data(), ld.data(), lg.data(), li.data(), &b.mat(), k);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j <= i; ++j)
                {
                    const double av = blasfeo_matel_wrap(&a.mat(), i, j);
                    const double bv = blasfeo_matel_wrap(&b.mat(), i, j);
                    const double error = std::abs(av - bv) / (1.0 + std::abs(av));
                    if (!std::isfinite(error)) return false;
                    worst = std::max(worst, error);
                }
            checker.check_stage(k, u, x, next, fd);
        }
    }
    std::cout << P::name << " checks=" << fd.checks << " failures=" << fd.failures
              << " reference_hessian_error=" << worst << '\n';
    return fd.failures == 0 && worst < 1e-10;
}
int main()
{
    bool ok = true;
#define CHECK(P) ok = check<P>() && ok;
    FATROP_BENCHMARK_PROBLEMS(CHECK)
#undef CHECK
    return ok ? 0 : 1;
}
