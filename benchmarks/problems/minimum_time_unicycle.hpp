#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
class MinimumTimeUnicycle : public ProblemBase<MinimumTimeUnicycle, 3, 3, 31>
{
public:
    static constexpr const char *name = "minimum-time-unicycle";
    static constexpr int terminal_equalities = 3;
    double target(int i) const { return i == 2 ? 1.5707963267948966 : 1.0; }
    double lower_control(int j) const { return j == 0 ? 0.0 : j == 1 ? -1.0 : 0.02; }
    double upper_control(int j) const { return j == 2 ? 0.2 : 1.0; }
    double initial_control(int, int j) const { return j == 2 ? 0.1 : 0.5; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    { return i == 0 ? x[0] + u[2] * u[0] * ad::cos(x[2])
             : i == 1 ? x[1] + u[2] * u[0] * ad::sin(x[2])
                      : x[2] + u[2] * u[1]; }
    template <class T> T cost(const T *u, const T *, int k) const
    { return minimum_time_cost(u, k, 0.02); }
};
}
