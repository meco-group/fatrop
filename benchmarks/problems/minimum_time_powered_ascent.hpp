#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
class MinimumTimePoweredAscent : public ProblemBase<MinimumTimePoweredAscent, 3, 2, 41>
{
public:
    static constexpr const char *name = "minimum-time-powered-ascent";
    static constexpr int terminal_equalities = 2;
    double initial(int i) const { return i == 2 ? 1.0 : 0.0; }
    double target(int i) const { return i == 0 ? 1.0 : 0.0; }
    double lower_control(int j) const { return j == 1 ? 0.02 : 0.0; }
    double upper_control(int j) const { return j == 1 ? 0.2 : 2.0; }
    double initial_control(int, int j) const { return j == 1 ? 0.1 : 0.0; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    { return i == 0 ? x[0] + u[1] * x[1]
             : i == 1 ? x[1] + u[1] * (u[0] * ad::exp(-x[2]) - 0.4
                            - 0.05 * x[1] * x[1] * x[1])
                      : x[2] - 0.015 * u[1] * u[0] * u[0]; }
    template <class T> T cost(const T *u, const T *, int k) const
    { return minimum_time_cost(u, k); }
};
}
