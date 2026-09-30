#pragma once
#include "problem_base.hpp"
namespace fatrop_benchmark
{
template <class Derived> class GoddardModel : public ProblemBase<Derived, 3, 1, 81>
{
public:
    double initial(int i) const { return i == 2 ? 1.0 : 0.0; }
    double lower_control(int) const { return 0.0; }
    double upper_control(int) const { return 1.0; }
    double initial_control(int, int) const { return 0.6; }
    int extra_inequalities(int k) const
    { return (Derived::max_q ? 1 : 0) + (k + 1 == 81 ? 1 : 0); }
    double lower_inequality(int q, int) const
    { return Derived::max_q && q == 0 ? -std::numeric_limits<double>::infinity() : 0.5; }
    double upper_inequality(int q, int) const
    { return Derived::max_q && q == 0 ? 0.5 : std::numeric_limits<double>::infinity(); }
    template <class T> T inequality(int q, const T *, const T *x, int) const
    { return Derived::max_q && q == 0 ? ad::exp(-x[0]) * x[1] * x[1] : x[2]; }
    template <class T> T dynamics(int i, const T *u, const T *x, int) const
    { return i == 0 ? x[0] + 0.05 * x[1]
             : i == 1 ? x[1] + 0.05 * (u[0] / x[2] - 0.2
                           - 0.04 * ad::exp(-x[0]) * x[1] * x[1])
                      : x[2] - 0.01 * u[0]; }
    template <class T> T cost(const T *u, const T *x, int k) const
    {
        T result = this->get_nu(k) ? 0.001 * u[0] * u[0] : T(0.0);
        return k + 1 == 81 ? result - x[0] : result;
    }
};
}
