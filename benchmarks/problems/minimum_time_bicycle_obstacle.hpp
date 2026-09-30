#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
class MinimumTimeBicycleObstacle : public ProblemBase<MinimumTimeBicycleObstacle, 4, 3, 61>
{
public:
    static constexpr const char *name = "minimum-time-bicycle-obstacle";
    static constexpr int terminal_equalities = 4;
    double initial(int i) const { return i == 0 ? -1.0 : 0.0; }
    double target(int i) const { return i == 0 ? 1.0 : 0.0; }
    double lower_control(int j) const { return j == 2 ? 0.02 : j == 1 ? -0.5 : -1.0; }
    double upper_control(int j) const { return j == 2 ? 0.15 : j == 1 ? 0.5 : 1.0; }
    double initial_control(int, int j) const { return j == 2 ? 0.1 : 0.0; }
    double initial_state(int k, int i) const
    {
        constexpr double pi = 3.141592653589793;
        double s = double(k) / 60.0, duration = 6.0;
        double dx = 12.0 * s * (1.0 - s);
        double dy = 0.45 * pi * std::sin(2.0 * pi * s);
        if (i == 0) return -1.0 + 2.0 * (3.0 * s * s - 2.0 * s * s * s);
        if (i == 1) return 0.45 * std::sin(pi * s) * std::sin(pi * s);
        if (i == 2) return s == 0.0 || s == 1.0 ? 0.0 : std::atan2(dy, dx);
        return std::sqrt(dx * dx + dy * dy) / duration;
    }
    int extra_inequalities(int) const { return 1; }
    double lower_inequality(int, int) const { return 0.04; }
    template <class T> T inequality(int, const T *, const T *x, int) const
    { return x[0] * x[0] + x[1] * x[1]; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    { return i == 0 ? x[0] + u[2] * x[3] * ad::cos(x[2])
             : i == 1 ? x[1] + u[2] * x[3] * ad::sin(x[2])
             : i == 2 ? x[2] + u[2] * x[3] * ad::tan(u[1])
                      : x[3] + u[2] * u[0]; }
    template <class T> T cost(const T *u, const T *, int k) const
    { return minimum_time_cost(u, k); }
};
}
