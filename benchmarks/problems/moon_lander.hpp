#pragma once
#include "moon_lander_model.hpp"
namespace fatrop_benchmark
{
class MoonLander : public MoonLanderModel<MoonLander, 2, 81>
{
public:
    static constexpr const char *name = "moon-lander";
    static constexpr int terminal_equalities = 4;
};
}
