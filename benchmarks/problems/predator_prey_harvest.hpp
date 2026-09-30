#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
class PredatorPreyHarvest : public ProblemBase<PredatorPreyHarvest, 2, 1, 70>
{
public:
    static constexpr const char *name = "predator-prey-harvest";
    static constexpr int terminal_equalities = 0;
    double initial(int i) const { return i == 0 ? 1.2 : 0.6; }
    double target(int) const { return 0.8; }
    double lower_control(int) const { return 0.0; }
    double upper_control(int) const { return 1.0; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    { return i == 0 ? x[0] + 0.1 * (0.5 * x[0] - 0.6 * x[0] * x[1] - u[0] * x[0])
                    : x[1] + 0.1 * (-0.4 * x[1] + 0.5 * x[0] * x[1]); }
    template <class T> T cost(const T *u, const T *x, int k) const
    { return quadratic_cost(u, x, k); }
};
}
