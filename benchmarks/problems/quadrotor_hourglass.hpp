#pragma once
#include "quadrotor_model.hpp"
namespace fatrop_benchmark
{
class QuadrotorHourglass : public QuadrotorModel<QuadrotorHourglass, 2>
{
public:
    static constexpr const char *name = "quadrotor-hourglass";
    static constexpr int terminal_equalities = 6;
    static constexpr bool high_arc = false;
    int extra_inequalities(int) const { return 1; }
    double upper_inequality(int, int) const { return 0.0; }
    template <class T> T inequality(int, const T *, const T *x, int) const
    { return quadrotor_hourglass_value(x); }
};
}
