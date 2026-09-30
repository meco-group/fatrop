#pragma once
#include "tethered_payload.hpp"
namespace fatrop_benchmark
{
class MinimumTimeTetheredPayload : public ProblemBase<MinimumTimeTetheredPayload, 4, 2, 81>
{
public:
    static constexpr const char *name = "minimum-time-tethered-payload";
    static constexpr int terminal_equalities = 4;
    double initial(int i) const { return i == 2 ? 1.5 : 0.0; }
    double target(int i) const { return i == 0 ? 1.0 : i == 2 ? 2.0 : 0.0; }
    double lower_control(int j) const { return j == 1 ? 0.02 : -2.0; }
    double upper_control(int j) const { return j == 1 ? 0.1 : 2.0; }
    double initial_control(int, int j) const { return j == 1 ? 0.1 : 0.0; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    {
        T spring = tethered_spring(x);
        return i == 0 ? x[0] + u[1] * x[1]
             : i == 1 ? x[1] + u[1] * (u[0] + spring - 0.1 * x[1])
             : i == 2 ? x[2] + u[1] * x[3]
                      : x[3] + u[1] * (-spring - 0.1 * x[3]);
    }
    template <class T> T cost(const T *u, const T *, int k) const
    { return minimum_time_cost(u, k); }
};
}
