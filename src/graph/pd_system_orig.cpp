//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/pd_system_orig.hpp"

#include "fatrop/graph/hessian.hpp"
#include "fatrop/graph/jacobian.hpp"
#include "fatrop/graph/problem_info.hpp"
#include "fatrop/linear_algebra/blasfeo_operations.hpp"

using namespace fatrop;

LinearSystem<PdSystemType<GraphProblem>>::LinearSystem(
    const ProblemInfo<GraphProblem> &info, Jacobian<GraphProblem> &jac,
    Hessian<GraphProblem> &hess, const VecRealView &D_x, bool De_is_zero, const VecRealView &D_e,
    const VecRealView &Sl_i, const VecRealView &Su_i, const VecRealView &Zl_i,
    const VecRealView &Zu_i, VecRealView &rhs_f_x, VecRealView &rhs_f_s, VecRealView &rhs_g,
    VecRealView &rhs_cl, VecRealView &rhs_cu)
    : info_(info), m_(info.number_of_tangent_variables + 3 * info.number_of_slack_variables +
                      info.number_of_eq_constraints),
      jac_(jac), hess_(hess), D_x_(D_x), De_is_zero_(De_is_zero), D_e_(D_e), Sl_i_(Sl_i),
      Su_i_(Su_i), Zl_i_(Zl_i), Zu_i_(Zu_i), rhs_f_x_(rhs_f_x), rhs_f_s_(rhs_f_s), rhs_g_(rhs_g),
      rhs_cl_(rhs_cl), rhs_cu_(rhs_cu)
{
}

Index LinearSystem<PdSystemType<GraphProblem>>::m(const ProblemInfo<GraphProblem> &info)
{
    return info.number_of_tangent_variables + 3 * info.number_of_slack_variables +
           info.number_of_eq_constraints;
}

void LinearSystem<PdSystemType<GraphProblem>>::get_rhs(VecRealView &out)
{
    out.block(info_.number_of_tangent_variables, info_.pd_orig_offset_primal) = rhs_f_x_;
    out.block(info_.number_of_slack_variables, info_.pd_orig_offset_slack) = rhs_f_s_;
    out.block(info_.number_of_eq_constraints, info_.pd_orig_offset_mult) = rhs_g_;
    out.block(info_.number_of_slack_variables, info_.pd_orig_offset_zl) = rhs_cl_;
    out.block(info_.number_of_slack_variables, info_.pd_orig_offset_zu) = rhs_cu_;
}

void LinearSystem<PdSystemType<GraphProblem>>::set_rhs(const VecRealView &in)
{
    rhs_f_x_ = in.block(info_.number_of_tangent_variables, info_.pd_orig_offset_primal);
    rhs_f_s_ = in.block(info_.number_of_slack_variables, info_.pd_orig_offset_slack);
    rhs_g_ = in.block(info_.number_of_eq_constraints, info_.pd_orig_offset_mult);
    rhs_cl_ = in.block(info_.number_of_slack_variables, info_.pd_orig_offset_zl);
    rhs_cu_ = in.block(info_.number_of_slack_variables, info_.pd_orig_offset_zu);
}

void LinearSystem<PdSystemType<GraphProblem>>::apply_on_right(const VecRealView &x, Scalar alpha,
                                                              const VecRealView &y,
                                                              VecRealView &out)
{
    // Layout mirrors the dense PD system: see dense/pd_system_orig.cpp. Graph
    // problems have no path equality block — the equality multipliers vector
    // contains only the slack-equality block, so `mult_eq` has size 0.
    VecRealView x_primal = x.block(info_.number_of_tangent_variables, info_.pd_orig_offset_primal);
    VecRealView x_slack = x.block(info_.number_of_slack_variables, info_.pd_orig_offset_slack);
    VecRealView mult = x.block(info_.number_of_eq_constraints, info_.pd_orig_offset_mult);
    VecRealView mult_ineq = mult; // all multipliers are slack-eq (offset_g_eq_slack == 0)
    VecRealView zl = x.block(info_.number_of_slack_variables, info_.pd_orig_offset_zl);
    VecRealView zu = x.block(info_.number_of_slack_variables, info_.pd_orig_offset_zu);
    VecRealView out_x = out.block(info_.number_of_tangent_variables, info_.pd_orig_offset_primal);
    VecRealView out_s = out.block(info_.number_of_slack_variables, info_.pd_orig_offset_slack);
    VecRealView out_mult = out.block(info_.number_of_eq_constraints, info_.pd_orig_offset_mult);
    VecRealView out_mult_ineq = out_mult;
    VecRealView out_zl = out.block(info_.number_of_slack_variables, info_.pd_orig_offset_zl);
    VecRealView out_zu = out.block(info_.number_of_slack_variables, info_.pd_orig_offset_zu);
    VecRealView D_x_x = D_x_.block(info_.number_of_tangent_variables, 0);
    VecRealView D_x_s = D_x_.block(info_.number_of_slack_variables, info_.offset_slack);

    out = alpha * y;
    hess_.apply_on_right(info_, x_primal, 1.0, out_x, out_x);
    out_x = out_x + D_x_x * x_primal;
    jac_.transpose_apply_on_right(info_, mult, 1.0, out_x, out_x);
    out_s = out_s + D_x_s * x_slack;
    out_s = out_s - mult_ineq;
    out_s = out_s - zl + zu;
    jac_.apply_on_right(info_, x_primal, 1.0, out_mult, out_mult);
    out_mult_ineq = out_mult_ineq - x_slack;
    if (!De_is_zero_)
    {
        // Only slack-equality block carries a (non-zero) D_e contribution.
        out_mult_ineq = out_mult_ineq - D_e_ * mult_ineq;
    }
    out_zl = out_zl + Zl_i_ * x_slack;
    out_zl = out_zl + Sl_i_ * zl;
    out_zu = out_zu - Zu_i_ * x_slack;
    out_zu = out_zu + Su_i_ * zu;
}
