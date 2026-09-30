#pragma once
#include "problem_base.hpp"

namespace fatrop_benchmark
{
template <class Derived, int NU, int STAGES>
class MoonLanderModel : public ProblemBase<Derived, 5, NU, STAGES>
{
public:
    using Base = ProblemBase<Derived, 5, NU, STAGES>;
    double initial(int i) const
    { return i == 1 || i == 4 ? 1.0 : i == 2 ? 0.2 : i == 3 ? -0.1 : 0.0; }
    double target(int) const { return 0.0; }
    double lower_control(int j) const { return NU == 3 && j == 2 ? 0.01 : -1.0; }
    double upper_control(int j) const { return NU == 3 && j == 2 ? 0.1 : 1.0; }
    double initial_control(int, int j) const { return NU == 3 && j == 2 ? 0.1 : 0.0; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    {
        T step = NU == 3 ? u[2] : T(0.05);
        return i == 0 ? x[0] + step * x[2]
             : i == 1 ? x[1] + step * x[3]
             : i == 2 ? x[2] + step * u[0] / x[4]
             : i == 3 ? x[3] + step * (u[1] / x[4] - 0.15)
                      : x[4] - 0.03 * step * (u[0] * u[0] + u[1] * u[1]);
    }
    template <class T> T cost(const T *u, const T *x, int k) const
    {
        if constexpr (NU == 3) return Base::minimum_time_cost(u, k);
        T result = this->get_nu(k) ? 0.025 * (u[0] * u[0] + u[1] * u[1]) : T(0.0);
        return k + 1 == STAGES ? result - 0.2 * x[4] : result;
    }
};

template <class Derived, int NU, int STAGES>
class MoonObstacleModel : public MoonLanderModel<Derived, NU, STAGES>
{
public:
    int extra_inequalities(int) const { return 1; }
    double lower_inequality(int, int) const { return 0.04; }
    template <class T> T inequality(int, const T *, const T *x, int) const
    { T dy = x[1] - 0.5; return x[0] * x[0] + dy * dy; }
    double initial_state(int k, int i) const
    {
        double s = double(k) / double(STAGES - 1);
        constexpr double pi = 3.141592653589793;
        if (i == 0) return 0.35 * std::sin(pi * s);
        if (i == 1) return 1.0 - s;
        if (i == 2) return 0.35 * pi / 4.0 * std::cos(pi * s);
        if (i == 3) return -0.25;
        return 1.0;
    }
};
}
