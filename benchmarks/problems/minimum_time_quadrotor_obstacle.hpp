#pragma once
#include "quadrotor_model.hpp"
namespace fatrop_benchmark
{
class MinimumTimeQuadrotorObstacle : public QuadrotorModel<MinimumTimeQuadrotorObstacle, 3>
{
public:
    static constexpr const char *name = "minimum-time-quadrotor-obstacle";
    static constexpr int terminal_equalities = 6;
    static constexpr bool high_arc = true;
    int extra_inequalities(int) const { return 1; }
    double lower_inequality(int, int) const { return 0.09; }
    template <class T> T inequality(int, const T *, const T *x, int) const
    { return quadrotor_obstacle_value(x); }
};
}
