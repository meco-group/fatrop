#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
class CoupledChain : public ProblemBase<CoupledChain, 12, 3, 50>
{
public:
    static constexpr const char *name = "coupled-chain";
    static constexpr int terminal_equalities = 0;
    double target(int i) const { return i % 2 ? -0.5 : 1.0; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    {
        int previous = (i + states - 1) % states;
        T v = u[i % controls];
        return 0.82 * x[i] + 0.04 * x[previous] + 0.02 * x[i] * x[i]
             + 0.02 * x[i] * x[previous] + 0.18 * v + 0.015 * v * v;
    }
    template <class T> T cost(const T *u, const T *x, int k) const
    {
        T result = quadratic_cost(u, x, k);
        for (int i = 1; i < states; ++i)
        {
            T d = x[i] - x[i - 1];
            result = result + 0.04 * d * d;
        }
        return result;
    }
};
}
