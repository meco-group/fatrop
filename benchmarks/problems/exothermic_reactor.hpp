#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
class ExothermicReactor : public ProblemBase<ExothermicReactor, 2, 1, 81>
{
public:
    static constexpr const char *name = "exothermic-reactor";
    static constexpr int terminal_equalities = 0;
    double initial(int i) const { return i == 0 ? 1.2 : 0.4; }
    double target(int i) const { return i == 0 ? 0.3 : 0.6; }
    double lower_control(int) const { return -1.0; }
    double upper_control(int) const { return 1.0; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    {
        T rate = x[0] * x[0] * ad::exp(0.3 * x[1]);
        return i == 0 ? x[0] + 0.08 * (1.0 - x[0] - 0.5 * rate)
                      : x[1] + 0.08 * (-0.7 * x[1] + rate + u[0]);
    }
    template <class T> T cost(const T *u, const T *x, int k) const
    { return quadratic_cost(u, x, k); }
};
}
