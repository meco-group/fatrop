#pragma once
#include "goddard_model.hpp"
namespace fatrop_benchmark
{
class GoddardMaxQ : public GoddardModel<GoddardMaxQ>
{
public:
    static constexpr const char *name = "goddard-max-q";
    static constexpr int terminal_equalities = 0;
    static constexpr bool max_q = true;
};
}
