//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#include "fatrop/dense/pd_solver_resto.hpp"
#include "fatrop/dense/aug_system_solver.hpp"
#include "fatrop/dense/pd_solver_orig.hpp"
#include "fatrop/dense/pd_system_orig.hpp"
#include "fatrop/dense/pd_system_resto.hpp"
#include "fatrop/dense/problem_info.hpp"
#include "fatrop/linear_algebra/linear_solver.hxx"
#include "fatrop/linear_algebra/vector.hpp"

using namespace fatrop;

template class fatrop::LinearSolver<PdSolverResto<DenseType>, PdSystemResto<DenseType>>;

PdSolverResto<DenseType>::PdSolverResto(
    const ProblemInfo<DenseType> &info,
    const std::shared_ptr<PdSolverOrig<DenseType>> &orig_solver)
    : LinearSolver<PdSolverResto<DenseType>, PdSystemResto<DenseType>>(
          LinearSystem<PdSystemResto<DenseType>>::m(info)),
      orig_solver_(orig_solver), D_e_orig_(info.number_of_eq_constraints),
      rhs_g_orig_(info.number_of_eq_constraints), f_pp_(info.number_of_eq_constraints),
      f_nn_(info.number_of_eq_constraints), Xpm1_(info.number_of_eq_constraints),
      Xnm1_(info.number_of_eq_constraints),
      x_orig_(LinearSystem<PdSystemType<DenseType>>::m(info))
{
}

void PdSolverResto<DenseType>::reduce(LinearSystem<PdSystemResto<DenseType>> &ls)
{
    const ProblemInfo<DenseType> &info = ls.info_;
    VecRealView f_p = ls.rhs_f_s_.block(info.number_of_eq_constraints, info.offset_p);
    VecRealView f_n = ls.rhs_f_s_.block(info.number_of_eq_constraints, info.offset_n);
    VecRealView rhs_cp = ls.rhs_cl_.block(info.number_of_eq_constraints, info.offset_p);
    VecRealView rhs_cn = ls.rhs_cl_.block(info.number_of_eq_constraints, info.offset_n);
    VecRealView P = ls.Sl_i_.block(info.number_of_eq_constraints, info.offset_p);
    VecRealView N = ls.Sl_i_.block(info.number_of_eq_constraints, info.offset_n);
    VecRealView Zp = ls.Zl_i_.block(info.number_of_eq_constraints, info.offset_p);
    VecRealView Zn = ls.Zl_i_.block(info.number_of_eq_constraints, info.offset_n);
    VecRealView Dxp = ls.D_x_.block(info.number_of_eq_constraints, info.offset_slack_p);
    VecRealView Dxn = ls.D_x_.block(info.number_of_eq_constraints, info.offset_slack_n);

    f_pp_ = if_else(info.constraint_allows_dual_damping, f_p + 1.0 / P * rhs_cp,
                    VecRealScalar(f_pp_.m(), 0.));
    f_nn_ = if_else(info.constraint_allows_dual_damping, f_n + 1.0 / N * rhs_cn,
                    VecRealScalar(f_nn_.m(), 0.));

    Xpm1_ = if_else(info.constraint_allows_dual_damping, P / (P * Dxp + Zp),
                    VecRealScalar(Xpm1_.m(), 0.));
    Xnm1_ = if_else(info.constraint_allows_dual_damping, N / (N * Dxn + Zn),
                    VecRealScalar(Xnm1_.m(), 0.));

    D_e_orig_ = ls.D_e_ + Xpm1_ + Xnm1_;
    rhs_g_orig_ = ls.rhs_g_ - Xpm1_ * f_pp_ + Xnm1_ * f_nn_;
}

void PdSolverResto<DenseType>::dereduce(LinearSystem<PdSystemResto<DenseType>> &ls,
                                        VecRealView &x)
{
    const ProblemInfo<DenseType> &info = ls.info_;

    VecRealView orig_primal_x =
        x_orig_.block(info.number_of_tangent_variables, info.pd_orig_offset_primal);
    VecRealView orig_primal_s =
        x_orig_.block(info.number_of_slack_variables, info.pd_orig_offset_slack);
    VecRealView orig_mult = x_orig_.block(info.number_of_eq_constraints, info.pd_orig_offset_mult);
    VecRealView orig_zl = x_orig_.block(info.number_of_slack_variables, info.pd_orig_offset_zl);
    VecRealView orig_zu = x_orig_.block(info.number_of_slack_variables, info.pd_orig_offset_zu);
    x.block(info.number_of_tangent_variables, info.pd_resto_offset_primal) = orig_primal_x;
    x.block(info.number_of_slack_variables, info.pd_resto_offset_slack) = orig_primal_s;
    x.block(info.number_of_eq_constraints, info.pd_resto_offset_mult) = orig_mult;
    x.block(info.number_of_slack_variables, info.pd_resto_offset_zl) = orig_zl;
    x.block(info.number_of_slack_variables, info.pd_resto_offset_zu) = orig_zu;

    VecRealView x_s = x.block(info.number_of_slack_variables_resto, info.pd_resto_offset_slack);
    VecRealView x_p = x_s.block(info.number_of_eq_constraints, info.offset_p);
    VecRealView x_n = x_s.block(info.number_of_eq_constraints, info.offset_n);
    x_p = if_else(info.constraint_allows_dual_damping, Xpm1_ * (-f_pp_ - orig_mult),
                  VecRealScalar(x_p.m(), 0.));
    x_n = if_else(info.constraint_allows_dual_damping, Xnm1_ * (-f_nn_ + orig_mult),
                  VecRealScalar(x_n.m(), 0.));

    VecRealView P = ls.Sl_i_.block(info.number_of_eq_constraints, info.offset_p);
    VecRealView N = ls.Sl_i_.block(info.number_of_eq_constraints, info.offset_n);
    VecRealView rhs_cp = ls.rhs_cl_.block(info.number_of_eq_constraints, info.offset_p);
    VecRealView rhs_cn = ls.rhs_cl_.block(info.number_of_eq_constraints, info.offset_n);
    VecRealView Zp = ls.Zl_i_.block(info.number_of_eq_constraints, info.offset_p);
    VecRealView Zn = ls.Zl_i_.block(info.number_of_eq_constraints, info.offset_n);

    VecRealView x_zl = x.block(info.number_of_slack_variables_resto, info.pd_resto_offset_zl);
    VecRealView x_zp = x_zl.block(info.number_of_eq_constraints, info.offset_p);
    VecRealView x_zn = x_zl.block(info.number_of_eq_constraints, info.offset_n);

    x_zp = if_else(info.constraint_allows_dual_damping, 1.0 / P * (-rhs_cp - Zp * x_p),
                   VecRealScalar(x_zp.m(), 0.));
    x_zn = if_else(info.constraint_allows_dual_damping, 1.0 / N * (-rhs_cn - Zn * x_n),
                   VecRealScalar(x_zn.m(), 0.));
}

LinsolReturnFlag PdSolverResto<DenseType>::solve_once_impl(
    LinearSystem<PdSystemResto<DenseType>> &ls, VecRealView &x)
{
    const ProblemInfo<DenseType> &info = ls.info_;
    reduce(ls);
    VecRealView rhs_f_s = ls.rhs_f_s_.block(info.number_of_slack_variables, 0);
    VecRealView rhs_cl = ls.rhs_cl_.block(info.number_of_slack_variables, 0);
    VecRealView rhs_cu = ls.rhs_cu_.block(info.number_of_slack_variables, 0);

    LinearSystem<PdSystemType<DenseType>> ls_orig(
        ls.info_, ls.jac_, ls.hess_,
        ls.D_x_.block(info.number_of_tangent_variables + info.number_of_slack_variables, 0), false,
        D_e_orig_, ls.Sl_i_.block(info.number_of_slack_variables, 0),
        ls.Su_i_.block(info.number_of_slack_variables, 0),
        ls.Zl_i_.block(info.number_of_slack_variables, 0),
        ls.Zu_i_.block(info.number_of_slack_variables, 0), ls.rhs_f_x_, rhs_f_s, rhs_g_orig_,
        rhs_cl, rhs_cu);
    LinsolReturnFlag ret = orig_solver_->solve_once(ls_orig, x_orig_);
    dereduce(ls, x);
    return ret;
}

void PdSolverResto<DenseType>::solve_rhs_impl(LinearSystem<PdSystemResto<DenseType>> &ls,
                                              VecRealView &x)
{
    const ProblemInfo<DenseType> &info = ls.info_;
    reduce(ls);
    VecRealView rhs_f_s = ls.rhs_f_s_.block(info.number_of_slack_variables, 0);
    VecRealView rhs_cl = ls.rhs_cl_.block(info.number_of_slack_variables, 0);
    VecRealView rhs_cu = ls.rhs_cu_.block(info.number_of_slack_variables, 0);

    LinearSystem<PdSystemType<DenseType>> ls_orig(
        ls.info_, ls.jac_, ls.hess_,
        ls.D_x_.block(info.number_of_tangent_variables + info.number_of_slack_variables, 0), false,
        D_e_orig_, ls.Sl_i_.block(info.number_of_slack_variables, 0),
        ls.Su_i_.block(info.number_of_slack_variables, 0),
        ls.Zl_i_.block(info.number_of_slack_variables, 0),
        ls.Zu_i_.block(info.number_of_slack_variables, 0), ls.rhs_f_x_, rhs_f_s, rhs_g_orig_,
        rhs_cl, rhs_cu);
    orig_solver_->solve_rhs(ls_orig, x_orig_);
    dereduce(ls, x);
}
