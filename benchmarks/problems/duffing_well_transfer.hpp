#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
class DuffingWellTransfer : public ProblemBase<DuffingWellTransfer, 2, 1, 81>
{
public:
    static constexpr const char *name = "duffing-well-transfer";
    static constexpr int terminal_equalities = 2;
    double initial(int i) const { return i == 0 ? -1.0 : 0.0; }
    double target(int i) const { return i == 0 ? 1.0 : 0.0; }
    double state_weight(int k) const { return k + 1 == stages ? 2.0 : 0.2; }
    double control_weight() const { return 0.05; }
    double lower_control(int) const { return -0.6; }
    double upper_control(int) const { return 0.6; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    { return i == 0 ? x[0] + 0.1 * x[1]
                    : x[1] + 0.1 * (x[0] - x[0] * x[0] * x[0] - 0.2 * x[1] + u[0]); }
    template <class T> T cost(const T *u, const T *x, int k) const
    { return quadratic_cost(u, x, k); }
};
}
