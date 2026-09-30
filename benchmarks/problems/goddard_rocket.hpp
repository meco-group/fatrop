#pragma once
#include "goddard_model.hpp"
namespace fatrop_benchmark
{
class GoddardRocket : public GoddardModel<GoddardRocket>
{
public:
    static constexpr const char *name = "goddard-rocket";
    static constexpr int terminal_equalities = 0;
    static constexpr bool max_q = false;
};
}
