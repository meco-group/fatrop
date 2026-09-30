#pragma once
#include "moon_lander_model.hpp"
namespace fatrop_benchmark
{
class MinimumTimeMoonLanderObstacle : public MoonObstacleModel<MinimumTimeMoonLanderObstacle, 3, 41>
{
public:
    static constexpr const char *name = "minimum-time-moon-lander-obstacle";
    static constexpr int terminal_equalities = 4;
};
}
