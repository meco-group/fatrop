#pragma once
#include "moon_lander_model.hpp"
namespace fatrop_benchmark
{
class MoonLanderObstacle : public MoonObstacleModel<MoonLanderObstacle, 2, 81>
{
public:
    static constexpr const char *name = "moon-lander-obstacle";
    static constexpr int terminal_equalities = 4;
};
}
