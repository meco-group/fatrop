#pragma once
#include "problem_base.hpp"

namespace fatrop_benchmark
{
template <class T> T quadrotor_obstacle_value(const T *x)
{
    T dx = x[0] - 0.15, dz = x[1] - 0.32;
    return dx * dx + dz * dz;
}
template <class T> T quadrotor_hourglass_value(const T *x)
{
    T width = 0.12 + 0.5 * x[0] * x[0];
    return x[1] * x[1] - width * width;
}

template <class Derived, int NU>
class QuadrotorModel : public ProblemBase<Derived, 6, NU, 81>
{
public:
    using Base = ProblemBase<Derived, 6, NU, 81>;
    double initial(int i) const { return i == 0 ? -1.0 : 0.0; }
    double target(int i) const { return i == 0 ? 1.0 : 0.0; }
    double lower_control(int j) const { return NU == 3 && j == 2 ? 0.02 : 0.0; }
    double upper_control(int j) const { return NU == 3 && j == 2 ? 0.1 : 1.2; }
    double initial_control(int, int j) const { return NU == 3 && j == 2 ? 0.1 : 0.25; }
    double initial_state(int k, int i) const
    {
        constexpr double pi = 3.141592653589793;
        const double s = double(k) / 80.0, duration = 6.4;
        if (i == 0) return -1.0 + 2.0 * (3.0 * s * s - 2.0 * s * s * s);
        if (i == 1) return Derived::high_arc ? 0.8 * std::sin(pi * s) : 0.0;
        if (i == 2) return 12.0 * s * (1.0 - s) / duration;
        if (i == 3) return Derived::high_arc ? 0.8 * pi * std::cos(pi * s) / duration : 0.0;
        if (i == 4) return -0.5 * std::sin(2.0 * pi * s);
        return -pi * std::cos(2.0 * pi * s) / duration;
    }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    {
        T step = NU == 3 ? u[2] : T(0.08);
        T thrust = u[0] + u[1];
        return i == 0 ? x[0] + step * x[2]
             : i == 1 ? x[1] + step * x[3]
             : i == 2 ? x[2] - step * thrust * ad::sin(x[4])
             : i == 3 ? x[3] + step * (thrust * ad::cos(x[4]) - 0.5)
             : i == 4 ? x[4] + step * x[5]
                      : x[5] + 2.0 * step * (u[1] - u[0]);
    }
    template <class T> T cost(const T *u, const T *x, int k) const
    {
        if constexpr (NU == 3) return Base::minimum_time_cost(u, k);
        return Base::quadratic_cost(u, x, k);
    }
};
}
