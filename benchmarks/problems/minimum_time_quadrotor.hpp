#pragma once
#include "quadrotor_model.hpp"
namespace fatrop_benchmark
{
class MinimumTimeQuadrotor : public QuadrotorModel<MinimumTimeQuadrotor, 3>
{
public:
    static constexpr const char *name = "minimum-time-quadrotor";
    static constexpr int terminal_equalities = 6;
    static constexpr bool high_arc = false;
};
}
