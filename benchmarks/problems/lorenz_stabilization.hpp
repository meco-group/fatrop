#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
class LorenzStabilization : public ProblemBase<LorenzStabilization, 3, 1, 61>
{
public:
    static constexpr const char *name = "lorenz-stabilization";
    static constexpr int terminal_equalities = 0;
    double initial(int) const { return 0.2; }
    double lower_control(int) const { return -10.0; }
    double upper_control(int) const { return 10.0; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    { return i == 0 ? x[0] + 0.2 * (x[1] - x[0])
             : i == 1 ? x[1] + 0.02 * (x[0] * (28.0 - x[2]) - x[1] + u[0])
                      : x[2] + 0.02 * (x[0] * x[1] - (8.0 / 3.0) * x[2]); }
    template <class T> T cost(const T *u, const T *x, int k) const
    { return quadratic_cost(u, x, k); }
};
}
