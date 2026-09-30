#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
class NonlinearPathLimit : public ProblemBase<NonlinearPathLimit, 2, 1, 50>
{
public:
    static constexpr const char *name = "nonlinear-path-limit";
    static constexpr int terminal_equalities = 0;
    double target(int i) const { return i == 0 ? 1.0 : 0.4; }
    int extra_inequalities(int) const { return 1; }
    double upper_inequality(int, int) const { return 0.5; }
    template <class T> T inequality(int, const T *u, const T *x, int k) const
    { return x[0] * x[0] + (get_nu(k) ? u[0] * u[0] : T(0.0)); }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    { return i == 0 ? 0.85 * x[0] + 0.25 * u[0]
                    : 0.8 * x[1] + 0.15 * x[0] * x[0]; }
    template <class T> T cost(const T *u, const T *x, int k) const
    { return quadratic_cost(u, x, k); }
};
}
