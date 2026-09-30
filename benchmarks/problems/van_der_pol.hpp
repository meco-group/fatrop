#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
class VanDerPol : public ProblemBase<VanDerPol, 2, 1, 60>
{
public:
    static constexpr const char *name = "van-der-pol";
    static constexpr int terminal_equalities = 0;
    double initial(int i) const { return i == 0 ? 1.5 : 0.0; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    { return i == 0 ? x[0] + 0.08 * x[1]
                    : x[1] + 0.08 * (1.2 * (1.0 - x[0] * x[0]) * x[1] - x[0] + u[0]); }
    template <class T> T cost(const T *u, const T *x, int k) const
    { return quadratic_cost(u, x, k); }
};
}
