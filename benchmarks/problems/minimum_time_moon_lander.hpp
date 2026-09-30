#pragma once
#include "moon_lander_model.hpp"
namespace fatrop_benchmark
{
class MinimumTimeMoonLander : public MoonLanderModel<MinimumTimeMoonLander, 3, 41>
{
public:
    static constexpr const char *name = "minimum-time-moon-lander";
    static constexpr int terminal_equalities = 4;
};
}
