#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
class Pendulum : public ProblemBase<Pendulum, 2, 1, 60>
{
public:
    static constexpr const char *name = "pendulum";
    static constexpr int terminal_equalities = 0;
    double initial(int i) const { return i == 0 ? 2.0 : 0.0; }
    double lower_control(int) const { return -2.0; }
    double upper_control(int) const { return 2.0; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    { return i == 0 ? x[0] + 0.12 * x[1]
                    : x[1] + 0.12 * (-ad::sin(x[0]) - 0.15 * x[1] + u[0]); }
    template <class T> T cost(const T *u, const T *x, int k) const
    {
        T result = get_nu(k) ? 0.04 * u[0] * u[0] : T(0.0);
        return result + state_weight(k) * (1.0 - ad::cos(x[0]))
            + 0.1 * state_weight(k) * x[1] * x[1];
    }
};
}
