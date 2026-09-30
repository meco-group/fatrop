#pragma once

#include <fatrop/fatrop.hpp>
#include "adjoint_ad.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace fatrop_benchmark
{
template <class Derived, int NX, int NU, int STAGES> class ProblemBase : public fatrop::OcpAbstract
{
public:
    using Index = fatrop::Index;
    using Scalar = fatrop::Scalar;
    static constexpr int states = NX, controls = NU, stages = STAGES;
    static_assert(NX + NU <= ad::capacity, "increase AD capacity for this problem");

    const Derived &model() const { return static_cast<const Derived &>(*this); }
    Index get_nx(Index) const override { return NX; }
    Index get_nu(Index k) const override { return k + 1 < STAGES ? NU : 0; }
    Index get_ng(Index k) const override
    { return k == 0 ? NX : k + 1 == STAGES ? Derived::terminal_equalities : 0; }
    Index get_ng_ineq(Index k) const override
    { return get_nu(k) + model().extra_inequalities(k); }
    Index get_horizon_length() const override { return STAGES; }

    double initial(int) const { return 0.0; }
    double target(int) const { return 0.0; }
    double lower_control(int) const { return -2.5; }
    double upper_control(int) const { return 2.5; }
    double initial_state(int, int i) const { return model().initial(i); }
    double initial_control(int, int) const { return 0.0; }
    int extra_inequalities(int) const { return 0; }
    template <class T> T inequality(int, const T *, const T *, int) const { return T(0.0); }
    double lower_inequality(int, int) const { return -std::numeric_limits<double>::infinity(); }
    double upper_inequality(int, int) const { return std::numeric_limits<double>::infinity(); }
    double state_weight(int k) const { return k + 1 == STAGES ? 8.0 : 1.0; }
    double control_weight() const { return 0.12; }
    template <class T> T quadratic_cost(const T *u, const T *x, int k) const
    {
        T result(0.0);
        for (int j = 0; j < get_nu(k); ++j)
            result = result + 0.5 * model().control_weight() * u[j] * u[j];
        for (int i = 0; i < NX; ++i)
        {
            T error = x[i] - model().target(i);
            result = result + 0.5 * model().state_weight(k) * error * error;
        }
        return result;
    }
    template <class T> T minimum_time_cost(const T *u, int k, double weight = 0.01) const
    {
        if (!get_nu(k)) return T(0.0);
        T result = u[NU - 1];
        for (int j = 0; j + 1 < NU; ++j)
            result = result + 0.5 * weight * u[j] * u[j];
        return result;
    }

    template <class T> T stage_cost(const T *u, const T *x, int k) const
    { return model().template cost<T>(u, x, k); }
    double dynamics_value(int i, const double *u, const double *x, int k) const
    { return model().template dynamics<double>(i, u, x, k); }
    double inequality_value(int i, const double *u, const double *x, int k) const
    { return model().template inequality<double>(i, u, x, k); }

    Index eval_BAbt(const Scalar *next, const Scalar *u, const Scalar *x,
                    MAT *r, Index k) override
    {
        fatrop::blasfeo_gese_wrap(r->m, r->n, 0.0, r, 0, 0);
        const int nu = get_nu(k), n = nu + NX;
        using A = ad::Forward<NX + NU>;
        std::array<A, NU> au{};
        std::array<A, NX> ax{};
        for (int j = 0; j < NU; ++j)
            if (j < nu) au[j] = A::variable(u[j], j);
        for (int i = 0; i < NX; ++i) ax[i] = A::variable(x[i], nu + i);
        for (int i = 0; i < NX; ++i)
        {
            A value = model().template dynamics<A>(i, au.data(), ax.data(), k);
            for (int j = 0; j < n; ++j)
                fatrop::blasfeo_matel_wrap(r, j, i) = value.tangent[j];
            fatrop::blasfeo_matel_wrap(r, n, i) = value.value - next[i];
        }
        return 0;
    }
    Index eval_b(const Scalar *next, const Scalar *u, const Scalar *x,
                 Scalar *r, Index k) override
    {
        for (int i = 0; i < NX; ++i) r[i] = dynamics_value(i, u, x, k) - next[i];
        return 0;
    }
    Index eval_RSQrqt(const Scalar *scale, const Scalar *u, const Scalar *x,
                      const Scalar *lam_dyn, const Scalar *, const Scalar *lam_eq_ineq,
                      MAT *r, Index k) override
    {
        fatrop::blasfeo_gese_wrap(r->m, r->n, 0.0, r, 0, 0);
        const int nu = get_nu(k), n = nu + NX;
        using A = ad::HessianNumber;
        thread_local ad::HessianTape tape;
        tape.reset();
        std::array<A, NU> au{};
        std::array<A, NX> ax{};
        std::array<int, ad::capacity> variables{};
        for (int j = 0; j < nu; ++j)
        { au[j] = tape.variable(u[j]); variables[j] = au[j].node; }
        for (int i = 0; i < NX; ++i)
        { ax[i] = tape.variable(x[i]); variables[nu + i] = ax[i].node; }
        A lagrangian = *scale * model().template cost<A>(au.data(), ax.data(), k);
        if (k + 1 < STAGES)
            for (int i = 0; i < NX; ++i)
                lagrangian = lagrangian + lam_dyn[i] *
                    model().template dynamics<A>(i, au.data(), ax.data(), k);
        for (int j = 0; j < model().extra_inequalities(k); ++j)
            lagrangian = lagrangian + lam_eq_ineq[nu + j] *
                model().template inequality<A>(j, au.data(), ax.data(), k);
        const auto hessian = tape.hessian(lagrangian, variables.data(), n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j <= i; ++j)
                fatrop::blasfeo_matel_wrap(r, i, j) = hessian[i * ad::capacity + j];
        return 0;
    }
    Index eval_Ggt(const Scalar *, const Scalar *x, MAT *r, Index k) override
    {
        const int ng = get_ng(k);
        if (!ng) return 0;
        fatrop::blasfeo_gese_wrap(r->m, r->n, 0.0, r, 0, 0);
        const int nu = get_nu(k);
        for (int i = 0; i < ng; ++i)
        {
            fatrop::blasfeo_matel_wrap(r, nu + i, i) = 1.0;
            fatrop::blasfeo_matel_wrap(r, nu + NX, i) =
                x[i] - (k == 0 ? model().initial(i) : model().target(i));
        }
        return 0;
    }
    Index eval_g(const Scalar *, const Scalar *x, Scalar *r, Index k) override
    {
        for (int i = 0; i < get_ng(k); ++i)
            r[i] = x[i] - (k == 0 ? model().initial(i) : model().target(i));
        return 0;
    }
    Index eval_Ggt_ineq(const Scalar *u, const Scalar *x, MAT *r, Index k) override
    {
        if (!get_ng_ineq(k)) return 0;
        fatrop::blasfeo_gese_wrap(r->m, r->n, 0.0, r, 0, 0);
        const int nu = get_nu(k), n = nu + NX;
        for (int j = 0; j < nu; ++j)
        {
            fatrop::blasfeo_matel_wrap(r, j, j) = 1.0;
            fatrop::blasfeo_matel_wrap(r, n, j) = u[j];
        }
        using A = ad::Forward<NX + NU>;
        std::array<A, NU> au{};
        std::array<A, NX> ax{};
        for (int j = 0; j < NU; ++j)
            if (j < nu) au[j] = A::variable(u[j], j);
        for (int i = 0; i < NX; ++i) ax[i] = A::variable(x[i], nu + i);
        for (int q = 0; q < model().extra_inequalities(k); ++q)
        {
            A value = model().template inequality<A>(q, au.data(), ax.data(), k);
            for (int j = 0; j < n; ++j)
                fatrop::blasfeo_matel_wrap(r, j, nu + q) = value.tangent[j];
            fatrop::blasfeo_matel_wrap(r, n, nu + q) = value.value;
        }
        return 0;
    }
    Index eval_gineq(const Scalar *u, const Scalar *x, Scalar *r, Index k) override
    {
        const int nu = get_nu(k);
        for (int j = 0; j < nu; ++j) r[j] = u[j];
        for (int q = 0; q < model().extra_inequalities(k); ++q)
            r[nu + q] = inequality_value(q, u, x, k);
        return 0;
    }
    Index eval_rq(const Scalar *scale, const Scalar *u, const Scalar *x,
                  Scalar *r, Index k) override
    {
        const int nu = get_nu(k), n = nu + NX;
        using A = ad::Number<double>;
        thread_local ad::Tape<double> tape;
        tape.reset();
        std::array<A, NU> au{};
        std::array<A, NX> ax{};
        for (int j = 0; j < nu; ++j) au[j] = tape.variable(u[j]);
        for (int i = 0; i < NX; ++i) ax[i] = tape.variable(x[i]);
        A value = model().template cost<A>(au.data(), ax.data(), k);
        tape.backward(value);
        for (int j = 0; j < nu; ++j) r[j] = *scale * tape.adjoint(au[j]);
        for (int i = 0; i < NX; ++i) r[nu + i] = *scale * tape.adjoint(ax[i]);
        return 0;
    }
    Index eval_L(const Scalar *scale, const Scalar *u, const Scalar *x,
                 Scalar *r, Index k) override
    { *r = *scale * stage_cost(u, x, k); return 0; }
    Index get_bounds(Scalar *lower, Scalar *upper, Index k) const override
    {
        const int nu = get_nu(k);
        for (int j = 0; j < nu; ++j)
        { lower[j] = model().lower_control(j); upper[j] = model().upper_control(j); }
        for (int q = 0; q < model().extra_inequalities(k); ++q)
        {
            lower[nu + q] = model().lower_inequality(q, k);
            upper[nu + q] = model().upper_inequality(q, k);
        }
        return 0;
    }
    Index get_initial_xk(Scalar *x, Index k) const override
    { for (int i = 0; i < NX; ++i) x[i] = model().initial_state(k, i); return 0; }
    Index get_initial_uk(Scalar *u, Index k) const override
    { for (int j = 0; j < get_nu(k); ++j) u[j] = model().initial_control(k, j); return 0; }
};
}
