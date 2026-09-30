#pragma once
#include <array>
#include <chrono>
#include <cstdint>

namespace fatrop_benchmark
{
template <class Problem> class TimedProblem : public Problem
{
public:
    enum Callback { DynamicsJacobian, DynamicsValue, Hessian, EqualityJacobian,
                    EqualityValue, InequalityJacobian, InequalityValue,
                    CostGradient, CostValue, CallbackCount };
private:
    using Clock = std::chrono::steady_clock;
    struct Scope
    {
        TimedProblem &owner;
        Callback callback;
        Clock::time_point start = Clock::now();
        Scope(TimedProblem &owner, Callback callback) : owner(owner), callback(callback) {}
        ~Scope()
        {
            const auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
                Clock::now() - start).count();
            owner.evaluation_nanoseconds_ += duration;
            ++owner.evaluation_calls_;
            owner.callback_nanoseconds_[callback] += duration;
            ++owner.callback_calls_[callback];
        }
    };
    std::int64_t evaluation_nanoseconds_ = 0;
    std::int64_t evaluation_calls_ = 0;
    std::array<std::int64_t, CallbackCount> callback_nanoseconds_{};
    std::array<std::int64_t, CallbackCount> callback_calls_{};

public:
    void reset_timing()
    { evaluation_nanoseconds_ = 0; evaluation_calls_ = 0;
      callback_nanoseconds_.fill(0); callback_calls_.fill(0); }
    double evaluation_ms() const { return 1e-6 * double(evaluation_nanoseconds_); }
    std::int64_t evaluation_calls() const { return evaluation_calls_; }
    double callback_ms(Callback callback) const
    { return 1e-6 * double(callback_nanoseconds_[callback]); }
    std::int64_t callback_calls(Callback callback) const
    { return callback_calls_[callback]; }

    fatrop::Index eval_BAbt(const fatrop::Scalar *next, const fatrop::Scalar *u,
                            const fatrop::Scalar *x, MAT *r, fatrop::Index k) override
    { Scope timer(*this, DynamicsJacobian); return Problem::eval_BAbt(next, u, x, r, k); }
    fatrop::Index eval_b(const fatrop::Scalar *next, const fatrop::Scalar *u,
                         const fatrop::Scalar *x, fatrop::Scalar *r, fatrop::Index k) override
    { Scope timer(*this, DynamicsValue); return Problem::eval_b(next, u, x, r, k); }
    fatrop::Index eval_RSQrqt(const fatrop::Scalar *scale, const fatrop::Scalar *u,
                              const fatrop::Scalar *x, const fatrop::Scalar *lam_dyn,
                              const fatrop::Scalar *lam_eq, const fatrop::Scalar *lam_eq_ineq,
                              MAT *r, fatrop::Index k) override
    { Scope timer(*this, Hessian); return Problem::eval_RSQrqt(scale, u, x, lam_dyn, lam_eq, lam_eq_ineq, r, k); }
    fatrop::Index eval_Ggt(const fatrop::Scalar *u, const fatrop::Scalar *x,
                           MAT *r, fatrop::Index k) override
    { Scope timer(*this, EqualityJacobian); return Problem::eval_Ggt(u, x, r, k); }
    fatrop::Index eval_g(const fatrop::Scalar *u, const fatrop::Scalar *x,
                         fatrop::Scalar *r, fatrop::Index k) override
    { Scope timer(*this, EqualityValue); return Problem::eval_g(u, x, r, k); }
    fatrop::Index eval_Ggt_ineq(const fatrop::Scalar *u, const fatrop::Scalar *x,
                                MAT *r, fatrop::Index k) override
    { Scope timer(*this, InequalityJacobian); return Problem::eval_Ggt_ineq(u, x, r, k); }
    fatrop::Index eval_gineq(const fatrop::Scalar *u, const fatrop::Scalar *x,
                             fatrop::Scalar *r, fatrop::Index k) override
    { Scope timer(*this, InequalityValue); return Problem::eval_gineq(u, x, r, k); }
    fatrop::Index eval_rq(const fatrop::Scalar *scale, const fatrop::Scalar *u,
                          const fatrop::Scalar *x, fatrop::Scalar *r, fatrop::Index k) override
    { Scope timer(*this, CostGradient); return Problem::eval_rq(scale, u, x, r, k); }
    fatrop::Index eval_L(const fatrop::Scalar *scale, const fatrop::Scalar *u,
                         const fatrop::Scalar *x, fatrop::Scalar *r, fatrop::Index k) override
    { Scope timer(*this, CostValue); return Problem::eval_L(scale, u, x, r, k); }
};
}
