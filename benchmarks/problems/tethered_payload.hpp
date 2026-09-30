#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
template <class T> T tethered_spring(const T *x)
{ T d = x[2] - x[0] - 1.0; return 0.8 * d + 0.5 * d * d * d; }

class TetheredPayload : public ProblemBase<TetheredPayload, 4, 1, 81>
{
public:
    static constexpr const char *name = "tethered-payload";
    static constexpr int terminal_equalities = 4;
    double initial(int i) const { return i == 2 ? 1.5 : 0.0; }
    double target(int i) const { return i == 0 ? 1.0 : i == 2 ? 2.0 : 0.0; }
    double lower_control(int) const { return -2.0; }
    double upper_control(int) const { return 2.0; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    {
        T spring = tethered_spring(x);
        return i == 0 ? x[0] + 0.08 * x[1]
             : i == 1 ? x[1] + 0.08 * (u[0] + spring - 0.1 * x[1])
             : i == 2 ? x[2] + 0.08 * x[3]
                      : x[3] + 0.08 * (-spring - 0.1 * x[3]);
    }
    template <class T> T cost(const T *u, const T *x, int k) const
    { return quadratic_cost(u, x, k); }
};
}
