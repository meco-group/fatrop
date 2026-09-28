#include "fatrop/dense/dims.hpp"
#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/aug_system_solver.hpp"
#include "fatrop/dense/pd_solver_orig.hpp"
#include "fatrop/ip_algorithm/limited_memory_bfgs.hpp"
#include "fatrop/ip_algorithm/full_bfgs.hpp"
#include "fatrop/dense/bfgs_matrix_ops.hpp"
#include "fatrop/graph/bfgs_matrix_ops.hpp"
#include "fatrop/graph/dims.hpp"
#include "fatrop/graph/problem_info.hpp"
#include "fatrop/graph/jacobian.hpp"
#include "fatrop/nlp/low_rank_aug_system_solver.hpp"
#include <gtest/gtest.h>
#include <deque>

using namespace fatrop;

namespace
{
using DenseMatrix = std::vector<std::vector<Scalar>>;
struct Pair { std::vector<Scalar> s, y; };

DenseMatrix explicit_bfgs(Index n, Scalar sigma, const std::deque<Pair> &pairs)
{
    DenseMatrix B(n, std::vector<Scalar>(n, 0.));
    for (Index i = 0; i < n; ++i) B[i][i] = sigma;
    for (const auto &p : pairs)
    {
        std::vector<Scalar> bs(n, 0.);
        Scalar sy = 0., sbs = 0.;
        for (Index i = 0; i < n; ++i)
        {
            for (Index j = 0; j < n; ++j) bs[i] += B[i][j] * p.s[j];
            sy += p.s[i] * p.y[i];
            sbs += p.s[i] * bs[i];
        }
        for (Index i = 0; i < n; ++i)
            for (Index j = 0; j < n; ++j)
                B[i][j] += p.y[i] * p.y[j] / sy - bs[i] * bs[j] / sbs;
    }
    return B;
}

void compare(const LowRankHessian &model, const DenseMatrix &B)
{
    for (Index i = 0; i < model.W.m(); ++i)
        for (Index j = 0; j < model.W.m(); ++j)
        {
            Scalar a = i == j ? model.sigma : 0.;
            for (Index k = 0; k < model.rank(); ++k)
                a += LowRankHessian::sign(k) * model.W(i, k) * model.W(j, k);
            EXPECT_NEAR(a, B[i][j], 3e-11);
        }
}
}

TEST(LimitedMemoryBfgs, CurrentMultipliersFrozenPairsAndMemoryEviction)
{
    // f = sum (i+1)x_i^2/2; c = sum x_i^2/2; d = sum (i+2)x_i^2/2.
    // Changing multipliers tests the Jacobian snapshot and frozen historical y's.
    const Index n = 4;
    ProblemDims<DenseType> dims(n, 1, 1);
    ProblemInfo<DenseType> info(dims);
    Jacobian<DenseType> jac(dims);
    LimitedMemoryBfgs<DenseType> history(info);
    history.set_history(3);
    VecRealAllocated x(n), g(n), lambda(2);
    std::vector<Scalar> previous(n, 0.);
    std::deque<Pair> pairs;
    std::shared_ptr<const LowRankHessian> saved;
    DenseMatrix saved_dense;
    for (Index step = 0; step < 10; ++step)
    {
        lambda(0) = .1 + .13 * step;
        lambda(1) = .3 - .01 * step;
        jac.Gg_eqt = 0.; jac.Gg_ineqt = 0.;
        Pair pair{std::vector<Scalar>(n), std::vector<Scalar>(n)};
        Scalar sy = 0., ss = 0.;
        for (Index i = 0; i < n; ++i)
        {
            x(i) = .4 * std::sin(.7 * (step + 1) * (i + 1));
            g(i) = (i + 1.) * x(i);
            jac.Gg_eqt(i, 0) = x(i);
            jac.Gg_ineqt(i, 0) = (i + 2.) * x(i);
            pair.s[i] = x(i) - previous[i];
            pair.y[i] = (i + 1. + lambda(0) + (i + 2.) * lambda(1)) * pair.s[i];
            sy += pair.s[i] * pair.y[i];
            ss += pair.s[i] * pair.s[i];
            previous[i] = x(i);
        }
        if (step)
        {
            if (pairs.size() == 3) pairs.pop_front();
            pairs.push_back(pair);
        }
        auto model = history.update(info, x, g, jac, lambda);
        auto B = explicit_bfgs(n, step ? sy / ss : 1., pairs);
        compare(*model, B);
        EXPECT_EQ(history.history_size(), std::min(step, 3));
        if (step == 2) { saved = model; saved_dense = B; }
        if (step > 2) compare(*saved, saved_dense);
    }
}

TEST(FullBfgs, DenseMatrixSatisfiesSecantEquation)
{
    ProblemDims<DenseType> dims(3, 0, 0);
    ProblemInfo<DenseType> info(dims);
    Jacobian<DenseType> jac(dims);
    FullBfgs<DenseType> bfgs(info);
    Hessian<DenseType> hess(dims);
    VecRealAllocated x(3), g(3), lambda(0);
    x = 0.; g = 0.;
    bfgs.update(info, x, g, jac, lambda);
    x(0) = 1.; x(1) = 2.; x(2) = -1.;
    g(0) = 4.; g(1) = 7.; g(2) = -2.;
    bfgs.update(info, x, g, jac, lambda);
    bfgs.copy_to(hess);
    for (Index i = 0; i < 3; ++i)
    {
        Scalar product = 0.;
        for (Index j = 0; j < 3; ++j)
        {
            EXPECT_NEAR(hess.Hht(i, j), hess.Hht(j, i), 1e-12);
            product += hess.Hht(i, j) * x(j);
        }
        EXPECT_NEAR(product, g(i), 1e-12);
    }
}

TEST(FullBfgs, GraphUpdatesStoredEdgesOnly)
{
    ProblemDims<fatrop::GraphProblem> dims({2, 1, 2}, {{0, 1}, {1, 2}}, {0, 0, 0});
    ProblemInfo<fatrop::GraphProblem> info(dims);
    Jacobian<fatrop::GraphProblem> jac(dims);
    FullBfgs<fatrop::GraphProblem> bfgs(info);
    Hessian<fatrop::GraphProblem> hess(dims);
    VecRealAllocated x(5), g(5), lambda(0);
    x = 0.; g = 0.;
    bfgs.update(info, x, g, jac, lambda);
    x = 1.;
    for (Index i = 0; i < 5; ++i) g(i) = i + 2.;
    bfgs.update(info, x, g, jac, lambda);
    bfgs.copy_to(hess);
    EXPECT_TRUE(hess.H.has_block(1, 0));
    EXPECT_TRUE(hess.H.has_block(2, 1));
    EXPECT_FALSE(hess.H.has_block(2, 0));
    EXPECT_NEAR(hess.H.block(1, 0)(0, 0), .2, 1e-12);
    EXPECT_NEAR(hess.H.block(1, 0)(0, 1), .4, 1e-12);
    EXPECT_NEAR(hess.H.block(2, 1)(0, 0), .8, 1e-12);
    EXPECT_NEAR(hess.H.block(2, 1)(1, 0), 1., 1e-12);
}

TEST(LimitedMemoryBfgs, SkipsResetOnFollowingUpdateAndScaling)
{
    ProblemDims<DenseType> dims(1, 0, 0);
    ProblemInfo<DenseType> info(dims);
    Jacobian<DenseType> jac(dims);
    LimitedMemoryBfgs<DenseType> history(info);
    VecRealAllocated x(1), g(1), lambda(0);
    x = 0.; g = 0.;
    EXPECT_DOUBLE_EQ(history.update(info, x, g, jac, lambda)->sigma, 1.);
    x = 1.; g = 3.;
    EXPECT_DOUBLE_EQ(history.update(info, x, g, jac, lambda)->sigma, 3.);
    x = 2.; g = 2.;
    history.update(info, x, g, jac, lambda);
    EXPECT_EQ(history.history_size(), 1);
    EXPECT_EQ(history.skipped_updates(), 1);
    x = 3.; g = 1.;
    history.update(info, x, g, jac, lambda);
    EXPECT_EQ(history.history_size(), 1);
    EXPECT_EQ(history.skipped_updates(), 2);
    x = 4.; g = 8.;
    auto reset = history.update(info, x, g, jac, lambda);
    EXPECT_EQ(history.history_size(), 0);
    EXPECT_EQ(reset->rank(), 0);
    EXPECT_DOUBLE_EQ(reset->sigma, 1.);
    history.set_init_max(2.);
    x = 5.; g = 11.;
    EXPECT_DOUBLE_EQ(history.update(info, x, g, jac, lambda)->sigma, 2.);
    history.set_initialization("constant");
    history.set_init_val(.5);
    x = 6.; g = 14.;
    EXPECT_DOUBLE_EQ(history.update(info, x, g, jac, lambda)->sigma, .5);
    history.reset();
    auto initial = history.update(info, x, g, jac, lambda);
    EXPECT_EQ(initial->rank(), 0);
    EXPECT_DOUBLE_EQ(initial->sigma, .5);
}

TEST(LimitedMemoryBfgs, RejectsInvalidOptionsAndNonEuclideanDimensions)
{
    ProblemDims<DenseType> dims(2, 0, 0);
    ProblemInfo<DenseType> info(dims);
    LimitedMemoryBfgs<DenseType> history(info);
    EXPECT_THROW(history.set_approximation("sr1"), std::invalid_argument);
    EXPECT_THROW(history.set_history(0), std::invalid_argument);
    EXPECT_THROW(history.set_history(7), std::invalid_argument);
    EXPECT_THROW(LimitedMemoryBfgs<DenseType>(info, 0), std::invalid_argument);
    EXPECT_THROW(history.set_max_skipping(0), std::invalid_argument);
    EXPECT_THROW(history.set_initialization("invalid"), std::invalid_argument);
    EXPECT_THROW(history.set_init_val(0.), std::invalid_argument);
    EXPECT_THROW(history.set_init_max(std::numeric_limits<Scalar>::infinity()), std::invalid_argument);
    Jacobian<DenseType> jac(dims);
    VecRealAllocated x(2), gradient(1), lambda(0);
    EXPECT_THROW(history.update(info, x, gradient, jac, lambda), std::invalid_argument);
}

TEST(LimitedMemoryBfgs, SkippedSnapshotAndScalar2Initialization)
{
    ProblemDims<DenseType> dims(2, 0, 0);
    ProblemInfo<DenseType> info(dims);
    Jacobian<DenseType> jac(dims);
    LimitedMemoryBfgs<DenseType> history(info);
    history.set_initialization("scalar2");
    VecRealAllocated x(2), g(2), lambda(0);
    x = 0.; g = 0.;
    history.update(info, x, g, jac, lambda);
    x(0) = 1.; x(1) = 0.; g(0) = -1.; g(1) = 0.;
    history.update(info, x, g, jac, lambda); // skip negative curvature, but save derivatives
    x(0) = 2.; x(1) = 1.; g(0) = 1.; g(1) = 4.;
    auto model = history.update(info, x, g, jac, lambda);
    // s=(1,1), y=(2,4), so scalar2=20/6, not the ratio from the initial snapshot.
    EXPECT_NEAR(model->sigma, 20. / 6., 1e-14);
    EXPECT_EQ(history.history_size(), 1);
    EXPECT_EQ(history.skipped_updates(), 0);
    history.update(info, x, g, jac, lambda); // identical x is a tiny step
    EXPECT_EQ(history.skipped_updates(), 1);
    EXPECT_EQ(history.history_size(), 1);
}

TEST(LimitedMemoryBfgs, WoodburyMatchesFullConstrainedPrimalDualSystemAndReusesFactors)
{
    const Index n = 5;
    ProblemDims<DenseType> dims(n, 2, 2);
    ProblemInfo<DenseType> info(dims);
    Jacobian<DenseType> jac(dims);
    HessianWithLowRank<DenseType> compact(dims);
    Hessian<DenseType> explicit_hessian(dims);
    LimitedMemoryBfgs<DenseType> history(info);
    VecRealAllocated x(n), grad(n), lambda(info.number_of_eq_constraints);
    lambda = 0.;
    jac.Gg_eqt = 0.; jac.Gg_ineqt = 0.;
    jac.Gg_eqt(0, 0) = 1.; jac.Gg_eqt(2, 0) = -.4;
    jac.Gg_eqt(1, 1) = 1.; jac.Gg_eqt(3, 1) = .7;
    jac.Gg_ineqt(2, 0) = 1.; jac.Gg_ineqt(4, 1) = 1.;
    std::deque<Pair> pairs;
    std::vector<Scalar> previous(n, 0.);
    std::shared_ptr<const LowRankHessian> model;
    Scalar sigma = 1.;
    for (Index step = 0; step < 5; ++step)
    {
        Pair p{std::vector<Scalar>(n), std::vector<Scalar>(n)};
        Scalar sy = 0., ss = 0.;
        for (Index i = 0; i < n; ++i)
        {
            x(i) = std::sin((step + 1.) * (i + 1.));
            grad(i) = (i + 1.) * x(i);
            p.s[i] = x(i) - previous[i]; p.y[i] = (i + 1.) * p.s[i];
            sy += p.s[i] * p.y[i]; ss += p.s[i] * p.s[i];
            previous[i] = x(i);
        }
        if (step) { pairs.push_back(p); sigma = sy / ss; }
        model = history.update(info, x, grad, jac, lambda);
    }
    compact.set_diagonal(model->sigma); compact.low_rank = model;
    explicit_hessian.set_zero();
    auto B = explicit_bfgs(n, sigma, pairs);
    for (Index i = 0; i < n; ++i)
        for (Index j = 0; j < n; ++j) explicit_hessian.Hht(i, j) = B[i][j];
    VecRealAllocated dx(n + 2), de(4), sl(2), su(2), zl(2), zu(2);
    VecRealAllocated fx(n), fs(2), fg(4), cl(2), cu(2);
    dx = .03; de = .02; sl = 1.1; su = 1.7; zl = .8; zu = .2;
    auto augmented = std::make_shared<AugSystemSolver<DenseType>>(info);
    auto augmented_ref = std::make_shared<AugSystemSolver<DenseType>>(info);
    PdSolverOrig<DenseType> solver(
        info, std::make_shared<LowRankAugSystemSolver<DenseType>>(info, augmented));
    PdSolverOrig<DenseType> reference(info, augmented_ref);
    solver.set_min_it_ref(1);
    for (bool regularized : {false, true})
    {
        LinearSystem<PdSystemType<DenseType>> ls(info, jac, compact, dx, !regularized, de,
                                               sl, su, zl, zu, fx, fs, fg, cl, cu);
        LinearSystem<PdSystemType<DenseType>> ref(info, jac, explicit_hessian, dx, !regularized, de,
                                                sl, su, zl, zu, fx, fs, fg, cl, cu);
        VecRealAllocated rhs(ls.m()), solution(ls.m()), expected(ls.m()), residual(ls.m());
        for (Index repeat = 0; repeat < 3; ++repeat)
        {
            for (Index i = 0; i < rhs.m(); ++i) rhs(i) = std::cos(i + .4 * repeat);
            ls.set_rhs(rhs);
            auto ret = repeat == 0 ? solver.solve_in_place(ls) : solver.solve_in_place_rhs(ls);
            EXPECT_TRUE(ret == LinsolReturnFlag::SUCCESS || ret == LinsolReturnFlag::ITREF_INCREASE);
            ls.get_rhs(solution);
            ref.set_rhs(rhs);
            ASSERT_EQ(reference.solve_in_place(ref), LinsolReturnFlag::SUCCESS);
            ref.get_rhs(expected);
            for (Index i = 0; i < rhs.m(); ++i) EXPECT_NEAR(solution(i), expected(i), 1e-10);
            ref.set_rhs(rhs);
            ref.apply_on_right(solution, 1., rhs, residual);
            EXPECT_LT(norm_inf(residual), 1e-10);
        }
    }
}
