#pragma once
#include "problem_base.hpp"

namespace fatrop_benchmark
{
class RaceCarSpeedLimit : public ProblemBase<RaceCarSpeedLimit, 2, 1, 81>
{
public:
    static constexpr const char *name = "race-car-speed-limit";
    static constexpr int terminal_equalities = 1;
    double target(int i) const { return i == 0 ? 1.0 : 0.0; }
    double lower_control(int) const { return 0.0; }
    double upper_control(int) const { return 1.0; }
    double initial_control(int, int) const { return 0.5; }
    int extra_inequalities(int) const { return 1; }
    double upper_inequality(int, int) const { return 0.0; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    { return i == 0 ? x[0] + 0.05 * x[1] : x[1] + 0.05 * (u[0] - x[1]); }
    template <class T> T inequality(int, const T *, const T *x, int) const
    { return x[1] - 1.0 + 0.5 * ad::sin(6.283185307179586 * x[0]); }
    template <class T> T cost(const T *u, const T *, int k) const
    { return get_nu(k) ? 0.01 * u[0] * u[0] : T(0.0); }
};
}
