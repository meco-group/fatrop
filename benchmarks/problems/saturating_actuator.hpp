#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
class SaturatingActuator : public ProblemBase<SaturatingActuator, 6, 2, 50>
{
public:
    static constexpr const char *name = "saturating-actuator";
    static constexpr int terminal_equalities = 0;
    double target(int) const { return 0.8; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    { return 0.9 * x[i] + 0.03 * ad::sin(x[i]) + 0.25 * ad::tanh(u[i % controls]); }
    template <class T> T cost(const T *u, const T *x, int k) const
    { return quadratic_cost(u, x, k); }
};
}
