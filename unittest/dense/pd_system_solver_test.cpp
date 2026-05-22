//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//
#include "../random_matrix.hpp"
#include "fatrop/dense/aug_system_solver.hpp"
#include "fatrop/dense/dims.hpp"
#include "fatrop/dense/hessian.hpp"
#include "fatrop/dense/jacobian.hpp"
#include "fatrop/dense/pd_solver_orig.hpp"
#include "fatrop/dense/pd_system_orig.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/dense/type.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include <gtest/gtest.h>
#include <memory>

using namespace fatrop;

class DensePdTest : public ::testing::Test
{
protected:
    static constexpr Index nx_ = 12;
    static constexpr Index ng_ = 4;
    static constexpr Index ng_ineq_ = 5;

    ProblemDims<DenseType> dims{nx_, ng_, ng_ineq_};
    ProblemInfo<DenseType> info{dims};
    Jacobian<DenseType> jacobian{dims};
    Hessian<DenseType> hessian{dims};
    VecRealAllocated sl{info.number_of_slack_variables};
    VecRealAllocated su{info.number_of_slack_variables};
    VecRealAllocated zl{info.number_of_slack_variables};
    VecRealAllocated zu{info.number_of_slack_variables};
    VecRealAllocated rhs_cl{info.number_of_slack_variables};
    VecRealAllocated rhs_cu{info.number_of_slack_variables};
    VecRealAllocated rhs_x{info.number_of_primal_variables};
    VecRealAllocated rhs_g{info.number_of_eq_constraints};
    VecRealAllocated rhs_s{info.number_of_slack_variables};
    VecRealAllocated D_x{info.number_of_primal_variables + info.number_of_slack_variables};
    VecRealAllocated D_eq{info.number_of_eq_constraints};
    std::shared_ptr<AugSystemSolver<DenseType>> aug_solver{
        std::make_shared<AugSystemSolver<DenseType>>(info)};
    PdSolverOrig<DenseType> pd_solver{info, aug_solver};

    void SetUp() override
    {
        jacobian.Gg_eqt.block(nx_, ng_, 0, 0) = ::test::random_matrix(nx_, ng_);
        jacobian.Gg_ineqt.block(nx_, ng_ineq_, 0, 0) = ::test::random_matrix(nx_, ng_ineq_);
        hessian.Hht.block(nx_, nx_, 0, 0) = ::test::random_spd_matrix(nx_);

        for (Index i = 0; i < info.number_of_primal_variables; ++i)
            rhs_x(i) = 1.0 * i;
        for (Index i = 0; i < D_x.m(); ++i)
            D_x(i) = 0.1 + 0.05 * i;
        for (Index i = 0; i < info.number_of_eq_constraints; ++i)
        {
            rhs_g(i) = 1.0 * i;
            D_eq(i) = 1e-2 * (i + 1);
        }
        for (Index i = 0; i < info.number_of_slack_variables; ++i)
        {
            sl(i) = 1. + 0.1 * i;
            su(i) = 1. + 0.2 * i;
            zl(i) = 1. + 0.3 * i;
            zu(i) = 1. + 0.4 * i;
            rhs_cl(i) = 1. + 0.5 * i;
            rhs_cu(i) = 1. + 0.6 * i;
            rhs_s(i) = 1. + 0.7 * i;
        }
    }
};

TEST_F(DensePdTest, TestSolve)
{
    LinearSystem<PdSystemType<DenseType>> ls(info, jacobian, hessian, D_x, false, D_eq, sl, su,
                                             zl, zu, rhs_x, rhs_s, rhs_g, rhs_cl, rhs_cu);
    VecRealAllocated x_full(ls.m());
    VecRealAllocated rhs_save(ls.m());
    VecRealAllocated tmp(ls.m());
    ls.get_rhs(rhs_save);
    LinsolReturnFlag ret = pd_solver.solve_in_place(ls);
    EXPECT_EQ(ret, LinsolReturnFlag::SUCCESS);
    ls.get_rhs(x_full);
    ls.set_rhs(rhs_save);
    ls.apply_on_right(x_full, 1.0, rhs_save, tmp);
    EXPECT_NEAR(norm_inf(tmp), 0., 1e-6);
}

TEST_F(DensePdTest, TestSolveDegen)
{
    LinearSystem<PdSystemType<DenseType>> ls(info, jacobian, hessian, D_x, true, D_eq, sl, su,
                                             zl, zu, rhs_x, rhs_s, rhs_g, rhs_cl, rhs_cu);
    VecRealAllocated x_full(ls.m());
    VecRealAllocated rhs_save(ls.m());
    VecRealAllocated tmp(ls.m());
    ls.get_rhs(rhs_save);
    LinsolReturnFlag ret = pd_solver.solve_in_place(ls);
    EXPECT_EQ(ret, LinsolReturnFlag::SUCCESS);
    ls.get_rhs(x_full);
    ls.set_rhs(rhs_save);
    ls.apply_on_right(x_full, 1.0, rhs_save, tmp);
    EXPECT_NEAR(norm_inf(tmp), 0., 1e-6);
}
