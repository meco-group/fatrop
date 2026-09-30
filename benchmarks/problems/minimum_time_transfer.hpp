#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
class MinimumTimeTransfer : public ProblemBase<MinimumTimeTransfer, 2, 2, 31>
{
public:
    static constexpr const char *name = "minimum-time-transfer";
    static constexpr int terminal_equalities = 2;
    double target(int i) const { return i == 0 ? 1.0 : 0.0; }
    double lower_control(int j) const { return j == 1 ? 0.02 : -1.0; }
    double upper_control(int j) const { return j == 1 ? 0.2 : 1.0; }
    double initial_control(int, int j) const { return j == 1 ? 0.1 : 0.0; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    { return i == 0 ? x[0] + u[1] * x[1]
                    : x[1] + u[1] * (u[0] - 0.1 * x[1] * x[1] * x[1]); }
    template <class T> T cost(const T *u, const T *, int k) const
    { return minimum_time_cost(u, k, 0.001); }
};
}
