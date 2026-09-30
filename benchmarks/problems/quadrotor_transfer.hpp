#pragma once
#include "quadrotor_model.hpp"
namespace fatrop_benchmark
{
class QuadrotorTransfer : public QuadrotorModel<QuadrotorTransfer, 2>
{
public:
    static constexpr const char *name = "quadrotor-transfer";
    static constexpr int terminal_equalities = 6;
    static constexpr bool high_arc = false;
};
}
