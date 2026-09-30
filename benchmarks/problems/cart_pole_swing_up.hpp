#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
class CartPoleSwingUp : public ProblemBase<CartPoleSwingUp, 4, 1, 61>
{
public:
    static constexpr const char *name = "cart-pole-swing-up";
    static constexpr int terminal_equalities = 4;
    double initial(int i) const { return i == 2 ? 3.141592653589793 : 0.0; }
    double lower_control(int) const { return -2.0; }
    double upper_control(int) const { return 2.0; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    { return i == 0 ? x[0] + 0.1 * x[1]
             : i == 1 ? x[1] + 0.1 * (u[0] - 0.1 * x[1]
                              + 0.2 * ad::sin(x[2]) * x[3] * x[3])
             : i == 2 ? x[2] + 0.1 * x[3]
                      : x[3] + 0.1 * (ad::sin(x[2]) - u[0] * ad::cos(x[2]) - 0.1 * x[3]); }
    template <class T> T cost(const T *u, const T *x, int k) const
    {
        T result = get_nu(k) ? 0.06 * u[0] * u[0] : T(0.0);
        for (int i : {0, 1, 3}) result = result + 0.5 * state_weight(k) * x[i] * x[i];
        return result + state_weight(k) * (1.0 - ad::cos(x[2]));
    }
};
}
