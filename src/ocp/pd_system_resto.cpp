//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#include "fatrop/ocp/pd_system_resto.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/pd_system_orig.hpp"
#include "fatrop/ocp/problem_info.hpp"

using namespace fatrop;
LinearSystem<PdSystemResto<OcpType>>::LinearSystem(
    const ProblemInfo<OcpType> &info, Jacobian<OcpType> &jac, Hessian<OcpType> &hess,
    const VecRealView &D_x, bool De_is_zero, const VecRealView &D_e, const VecRealView &Sl_i,
    const VecRealView &Su_i, const VecRealView &Zl_i, const VecRealView &Zu_i, VecRealView &rhs_f_x,
    VecRealView &rhs_f_s, VecRealView &rhs_g, VecRealView &rhs_cl, VecRealView &rhs_cu)
    : info_(info), m_(LinearSystem<PdSystemResto<OcpType>>::m(info)), jac_(jac), hess_(hess),
      D_x_(D_x), De_is_zero_(De_is_zero), D_e_(D_e), Sl_i_(Sl_i), Su_i_(Su_i), Zl_i_(Zl_i),
      Zu_i_(Zu_i), rhs_f_x_(rhs_f_x), rhs_f_s_(rhs_f_s), rhs_g_(rhs_g), rhs_cl_(rhs_cl),
      rhs_cu_(rhs_cu)
{
}
Index LinearSystem<PdSystemResto<OcpType>>::m(const ProblemInfo<OcpType> &info)
{
    return info.number_of_primal_variables + info.number_of_eq_constraints +
           3 * info.number_of_slack_variables_resto /* grad_s rhs_zl rhs_zu*/;
}

void LinearSystem<PdSystemResto<OcpType>>::get_rhs(VecRealView &out)
{
    out.block(info_.number_of_primal_variables, info_.pd_resto_offset_primal) = rhs_f_x_;
    out.block(info_.number_of_slack_variables_resto, info_.pd_resto_offset_slack) = rhs_f_s_;
    out.block(info_.number_of_eq_constraints, info_.pd_resto_offset_mult) = rhs_g_;
    out.block(info_.number_of_slack_variables_resto, info_.pd_resto_offset_zl) = rhs_cl_;
    out.block(info_.number_of_slack_variables_resto, info_.pd_resto_offset_zu) = rhs_cu_;
}
void LinearSystem<PdSystemResto<OcpType>>::set_rhs(const VecRealView &in)
{
    rhs_f_x_ = in.block(info_.number_of_primal_variables, info_.pd_resto_offset_primal);
    rhs_f_s_ = in.block(info_.number_of_slack_variables_resto, info_.pd_resto_offset_slack);
    rhs_g_ = in.block(info_.number_of_eq_constraints, info_.pd_resto_offset_mult);
    rhs_cl_ = in.block(info_.number_of_slack_variables_resto, info_.pd_resto_offset_zl);
    rhs_cu_ = in.block(info_.number_of_slack_variables_resto, info_.pd_resto_offset_zu);
}

//    [ H + D_x    0        A_e^T  A_d^T  A_i^T    0     0  ] [ x   ] = [ -f_x ]
//    [    0     D_x          0      0  [-I+I-I]' -I     I  ] [ spn ] = [ -f_s ]
//    [ A_e      [+I -I]   -D_e      0     0       0     0  ] [ λ_e ] = [ -g_e ]
//    [ A_d   [0  +I -I]     0       0     0       0     0  ] [ λ_d ] = [ -g_d ]
//    [ A_i   [-I +I -I]     0       0     -D_e    0     0  ] [ λ_i ] = [ -g_i ]
//    [   0     Zl_i         0       0     0     Sl_i    0  ] [ zl  ] = [ -cl  ]
//    [   0    -Zu_i         0       0     0      0    Su_i ] [ zu  ] = [ -cu  ]
void LinearSystem<PdSystemResto<OcpType>>::apply_on_right(const VecRealView &x, Scalar alpha,
                                                          const VecRealView &y, VecRealView &out)
{
    VecRealView x_primal = x.block(info_.number_of_primal_variables, info_.pd_resto_offset_primal);
    VecRealView x_slack =
        x.block(info_.number_of_slack_variables_resto, info_.pd_resto_offset_slack);
    VecRealView mult = x.block(info_.number_of_eq_constraints, info_.pd_resto_offset_mult);
    VecRealView mult_eq = mult.block(info_.number_of_g_eq_path, info_.offset_g_eq_path);
    VecRealView mult_ineq = mult.block(info_.number_of_g_eq_slack, info_.offset_g_eq_slack);
    VecRealView mult_dyn = mult.block(info_.number_of_g_eq_dyn, info_.offset_g_eq_dyn);
    VecRealView zl = x.block(info_.number_of_slack_variables_resto, info_.pd_resto_offset_zl);
    VecRealView zu = x.block(info_.number_of_slack_variables_resto, info_.pd_resto_offset_zu);

    VecRealView out_x = out.block(info_.number_of_primal_variables, info_.pd_resto_offset_primal);
    VecRealView out_s =
        out.block(info_.number_of_slack_variables_resto, info_.pd_resto_offset_slack);
    VecRealView out_mult = out.block(info_.number_of_eq_constraints, info_.pd_resto_offset_mult);
    VecRealView out_mult_eq = out_mult.block(info_.number_of_g_eq_path, info_.offset_g_eq_path);
    VecRealView out_mult_ineq = out_mult.block(info_.number_of_g_eq_slack, info_.offset_g_eq_slack);
    VecRealView out_mult_dyn = out_mult.block(info_.number_of_g_eq_dyn, info_.offset_g_eq_dyn);
    VecRealView out_zl = out.block(info_.number_of_slack_variables_resto, info_.pd_resto_offset_zl);
    VecRealView out_zu = out.block(info_.number_of_slack_variables_resto, info_.pd_resto_offset_zu);

    VecRealView D_x_x = D_x_.block(info_.number_of_primal_variables, info_.offset_primal);
    VecRealView D_x_s = D_x_.block(info_.number_of_slack_variables_resto, info_.offset_slack);

    VecRealView s = x_slack.block(info_.number_of_slack_variables, info_.offset_s);
    VecRealView p = x_slack.block(info_.number_of_eq_constraints, info_.offset_p);
    VecRealView n = x_slack.block(info_.number_of_eq_constraints, info_.offset_n);

    VecRealView out_s_s = out_s.block(info_.number_of_slack_variables, info_.offset_s);
    VecRealView out_s_p = out_s.block(info_.number_of_eq_constraints, info_.offset_p);
    VecRealView out_s_n = out_s.block(info_.number_of_eq_constraints, info_.offset_n);

    out = alpha * y;
    // out_x += Hx
    hess_.apply_on_right(info_, x_primal, 1.0, out_x, out_x);
    // out_x += D_x @ x
    out_x = out_x + D_x_x * x_primal;
    // out_x += [A_e^T A_d^T A_i^T] @ mult
    jac_.transpose_apply_on_right(info_, mult, 1.0, out_x, out_x);
    // out_s += D_x @ s
    out_s = out_s + D_x_s * x_slack;
    // out_s -= I @ lam_I
    out_s_s = out_s_s - mult_ineq;

    out_s_p.block(info_.number_of_g_eq_path, info_.offset_g_eq_path) =
        out_s_p.block(info_.number_of_g_eq_path, info_.offset_g_eq_path) +
        mult.block(info_.number_of_g_eq_path, info_.offset_g_eq_path);
    out_s_n.block(info_.number_of_g_eq_path, info_.offset_g_eq_path) =
        out_s_n.block(info_.number_of_g_eq_path, info_.offset_g_eq_path) -
        mult.block(info_.number_of_g_eq_path, info_.offset_g_eq_path);

    out_s_p.block(info_.number_of_g_eq_slack, info_.offset_g_eq_slack) =
        out_s_p.block(info_.number_of_g_eq_slack, info_.offset_g_eq_slack) +
        mult.block(info_.number_of_g_eq_slack, info_.offset_g_eq_slack);
    out_s_n.block(info_.number_of_g_eq_slack, info_.offset_g_eq_slack) =
        out_s_n.block(info_.number_of_g_eq_slack, info_.offset_g_eq_slack) -
        mult.block(info_.number_of_g_eq_slack, info_.offset_g_eq_slack);

    // out_s -= -I @ zl + I @ zu
    out_s = out_s - zl + zu;
    // out_mult += [A_e; A_d; A_i] @ x
    jac_.apply_on_right(info_, x_primal, 1.0, out_mult, out_mult);
    // out_mult_ineq += -I @ s
    out_mult_ineq = out_mult_ineq - s;
    // out_mult_eq += -D_e @ lam_e
    if (!De_is_zero_)
        out_mult_eq =
            out_mult_eq - D_e_.block(info_.number_of_g_eq_path, info_.offset_g_eq_path) * mult_eq;
    // out_mult_ineq = -D_i @ lam_i
    out_mult_ineq =
        out_mult_ineq - D_e_.block(info_.number_of_g_eq_slack, info_.offset_g_eq_slack) * mult_ineq;
    // out_mult += p - n
    out_mult = out_mult + p - n;
    // out_zl += -Zl_i @ s
    out_zl = out_zl + Zl_i_ * x_slack;
    // out_zl += S_i @ zl
    out_zl = out_zl + Sl_i_ * zl;
    // out_zu += Zu_i @ zu
    out_zu = out_zu - Zu_i_ * x_slack;
    // out_zu += S_i @ zu
    out_zu = out_zu + Su_i_ * zu;
}