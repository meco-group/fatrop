//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

/**
 *
 * Todo this implementation is tightly coupled with the ocp problem type but actually the algorithm
 * implemented here is general for any problem type with a pd solver. It should be moved to the
 * ip_algorithm module in an hxx file. With a more general interface.
 *
 */

#include "fatrop/ocp/pd_solver_resto.hpp"
#include "fatrop/linear_algebra/linear_solver.hxx"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/ocp/aug_system_solver.hpp"
#include "fatrop/ocp/pd_solver_orig.hpp"
#include "fatrop/ocp/pd_system_orig.hpp"
#include "fatrop/ocp/pd_system_resto.hpp"
#include "fatrop/ocp/problem_info.hpp"
using namespace fatrop;

// instantiate the template class
template class fatrop::LinearSolver<PdSolverResto<OcpType>, PdSystemResto<OcpType>>;

PdSolverResto<OcpType>::PdSolverResto(const ProblemInfo<OcpType> &info,
                                      const std::shared_ptr<PdSolverOrig<OcpType>> &pd_solver_orig)
    : LinearSolver<PdSolverResto<OcpType>, PdSystemResto<OcpType>>(
          LinearSystem<PdSystemResto<OcpType>>::m(info)), orig_solver_(pd_solver_orig),
      D_e_orig_(info.number_of_eq_constraints), rhs_g_orig_(info.number_of_eq_constraints),
      f_pp_(info.number_of_eq_constraints), f_nn_(info.number_of_eq_constraints),
      Xpm1_(info.number_of_eq_constraints), Xnm1_(info.number_of_eq_constraints),
      x_orig_(LinearSystem<PdSystemType<OcpType>>::m(info))
{
}
//
// [ H    A^T |             ] [x]  = - [ f_x ]
// [ A   -D_e |  +I -I      ] [l]  = - [  g  ]
//     ------ | ---------
// [       I  | Dx     -I   ] [p]  = - [ f_p ]
// [      -I  |   Dx     -I ] [n]  = - [ f_n ]
// [          |  Zp    P    ] [zp] = - [ rhs_cp ]
// [          |     Zn    N ] [zn] = - [ rhs_cn ]

// Note that we omitted the slack variables and dual bounds, and associated equations of the
// original system. They dont affect the solution of this analysis.

// Eliminate zp and zn:
// zp = P^{-1} (-rhs_cp - Zp p)
// zn = N^{-1} (-rhs_cn - Zn n)

//
// [ H    A^T |         ] [x]  = - [ f_x  ]
// [ A   -D_e |  +I -I  ] [l]  = - [  g   ]
//     ------ | -------
// [       I  |  Xp     ] [p]  = - [ f_pp ]
// [      -I  |     Xn  ] [n]  = - [ f_nn ]

// here f_pp = f_p + P^{-1} rhs_cp
//      f_nn = f_n + N^{-1} rhs_cn
//      Xp = D_x + P^{-1} Zp
//      Xn = D_x + N^{-1} Zn

// eliminate p = Xp{-1} (-f_pp - l)
//           n = Xn{-1} (-f_nn + l)

//
// [ H    A^T  ] [x]  = - [ f_x  ]
// [ A   -D_ee ] [l]  = - [  gg  ]

// where D_ee = D_e + Xp^{-1} + Xn^{-1}
//       gg = g + Xp^{-1} f_pp + Xn^{-1} f_nn

//  X_p^{-1} = 1/ (D_x + P^{-1} Zp) = P / (P D_x + Zp)
//  X_n^{-1} = 1/ (D_x + N^{-1} Zn) = N / (N D_x + Zn)

// This is the same structure as the original system.
// Note that for OCP we dont add the p and n variables for the dynamics equations as they would
// change the structure of the system. While they are available we assume they are zero.
// This is achieved by setting the lower and upper bounds to infity in the OCP problem.
// Because of how fatrop threats these kinds of bounds the related dual variables will be one and
// the slack variables will always be zero. This ensures that it doesnt affect the solution of the
// system.

void PdSolverResto<OcpType>::reduce(LinearSystem<PdSystemResto<OcpType>> &ls)
{
    const ProblemInfo<OcpType> &info = ls.info_;
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


    f_pp_ = f_p + 1.0 / P * rhs_cp;
    f_nn_ = f_n + 1.0 / N * rhs_cn;
    VecRealView P_eq = P.block(info.number_of_g_eq_path, info.offset_g_eq_path);
    VecRealView N_eq = N.block(info.number_of_g_eq_path, info.offset_g_eq_path);
    VecRealView Zp_eq = Zp.block(info.number_of_g_eq_path, info.offset_g_eq_path);
    VecRealView Zn_eq = Zn.block(info.number_of_g_eq_path, info.offset_g_eq_path);
    VecRealView Dxp_eq = Dxp.block(info.number_of_g_eq_path, info.offset_g_eq_path);
    VecRealView Dxn_eq = Dxn.block(info.number_of_g_eq_path, info.offset_g_eq_path);
    VecRealView Xpm1_eq =Xpm1_.block(info.number_of_g_eq_path, info.offset_g_eq_path);
    VecRealView Xnm1_eq =Xnm1_.block(info.number_of_g_eq_path, info.offset_g_eq_path);

    Xpm1_eq = P_eq / (P_eq * Dxp_eq + Zp_eq);
    Xnm1_eq = N_eq / (N_eq * Dxn_eq + Zn_eq);

    VecRealView P_ineq = P.block(info.number_of_g_eq_slack, info.offset_g_eq_slack);
    VecRealView N_ineq = N.block(info.number_of_g_eq_slack, info.offset_g_eq_slack);
    VecRealView Zp_ineq = Zp.block(info.number_of_g_eq_slack, info.offset_g_eq_slack);
    VecRealView Zn_ineq = Zn.block(info.number_of_g_eq_slack, info.offset_g_eq_slack);
    VecRealView Dxp_ineq = Dxp.block(info.number_of_g_eq_slack, info.offset_g_eq_slack);
    VecRealView Dxn_ineq = Dxn.block(info.number_of_g_eq_slack, info.offset_g_eq_slack);
    VecRealView Xpm1_ineq =Xpm1_.block(info.number_of_g_eq_slack, info.offset_g_eq_slack);
    VecRealView Xnm1_ineq =Xnm1_.block(info.number_of_g_eq_slack, info.offset_g_eq_slack);

    Xpm1_ineq = P_ineq / (P_ineq * Dxp_ineq + Zp_ineq);
    Xnm1_ineq = N_ineq / (N_ineq * Dxn_ineq + Zn_ineq);

    D_e_orig_ = ls.D_e_ + Xpm1_ + Xnm1_;
    rhs_g_orig_ = ls.rhs_g_ + Xpm1_ * f_pp_ + Xnm1_ * f_nn_;
}
void PdSolverResto<OcpType>::dereduce(LinearSystem<PdSystemResto<OcpType>> &ls, VecRealView &x)
{
    const ProblemInfo<OcpType> &info = ls.info_;

    VecRealView orig_primal_x =
        x_orig_.block(info.number_of_primal_variables, info.pd_orig_offset_primal);
    VecRealView orig_primal_s = x_orig_.block(info.number_of_slack_variables, info.pd_orig_offset_slack);
    VecRealView orig_mult = x_orig_.block(info.number_of_eq_constraints, info.pd_orig_offset_mult);
    VecRealView orig_zl = x_orig_.block(info.number_of_slack_variables, info.pd_orig_offset_zl);
    VecRealView orig_zu = x_orig_.block(info.number_of_slack_variables, info.pd_orig_offset_zu);
    // set x
    x.block(info.number_of_primal_variables, info.pd_resto_offset_primal) = orig_primal_x;
    x.block(info.number_of_slack_variables, info.pd_resto_offset_slack) = orig_primal_s;
    x.block(info.number_of_eq_constraints, info.pd_resto_offset_mult) = orig_mult;
    x.block(info.number_of_slack_variables, info.pd_resto_offset_zl) = orig_zl;
    x.block(info.number_of_slack_variables, info.pd_resto_offset_zu) = orig_zu;

    VecRealView x_s =
        x.block(info.number_of_slack_variables_resto, info.pd_resto_offset_slack);
    VecRealView x_p = x_s.block(info.number_of_eq_constraints, info.offset_p);
    VecRealView x_n = x_s.block(info.number_of_eq_constraints, info.offset_n);
    x_p = Xpm1_ * (-f_pp_ - orig_mult);
    x_n = Xnm1_ * (-f_nn_ + orig_mult);

    VecRealView P = ls.Sl_i_.block(info.number_of_eq_constraints, info.offset_p);
    VecRealView N = ls.Sl_i_.block(info.number_of_eq_constraints, info.offset_n);
    VecRealView rhs_cp = ls.rhs_cl_.block(info.number_of_eq_constraints, info.offset_p);
    VecRealView rhs_cn = ls.rhs_cl_.block(info.number_of_eq_constraints, info.offset_n);
    VecRealView Zp = ls.Zl_i_.block(info.number_of_eq_constraints, info.offset_p);
    VecRealView Zn = ls.Zl_i_.block(info.number_of_eq_constraints, info.offset_n);

    VecRealView x_zl = x.block(info.number_of_slack_variables_resto, info.pd_resto_offset_zl);
    VecRealView x_zp = x_zl.block(info.number_of_eq_constraints, info.offset_p);
    VecRealView x_zn = x_zl.block(info.number_of_eq_constraints, info.offset_n);

    x_zp = 1.0 / P * (-rhs_cp - Zp * x_p);
    x_zn = 1.0 / N * (-rhs_cn - Zn * x_n);
}
LinsolReturnFlag PdSolverResto<OcpType>::solve_once_impl(LinearSystem<PdSystemResto<OcpType>> &ls,
                                                         VecRealView &x)
{
    const ProblemInfo<OcpType> &info = ls.info_;
    reduce(ls);
    LinsolReturnFlag ret;
    // set up the original linear system
    VecRealView rhs_f_s = ls.rhs_f_s_.block(info.number_of_slack_variables, 0);
    VecRealView rhs_cl = ls.rhs_cl_.block(info.number_of_slack_variables, 0);
    VecRealView rhs_cu = ls.rhs_cu_.block(info.number_of_slack_variables, 0);

    LinearSystem<PdSystemType<OcpType>> ls_orig(
        ls.info_, ls.jac_, ls.hess_,
        ls.D_x_.block(info.number_of_primal_variables + info.number_of_slack_variables, 0), false,
        D_e_orig_, ls.Sl_i_.block(info.number_of_slack_variables, 0),
        ls.Su_i_.block(info.number_of_slack_variables, 0),
        ls.Zl_i_.block(info.number_of_slack_variables, 0),
        ls.Zu_i_.block(info.number_of_slack_variables, 0), ls.rhs_f_x_, rhs_f_s, rhs_g_orig_,
        rhs_cl, rhs_cu);
    // solve ls orig and save in x_orig_
    ret = orig_solver_->solve_once(ls_orig, x_orig_);
    dereduce(ls, x);
    return ret;
}
void PdSolverResto<OcpType>::solve_rhs_impl(LinearSystem<PdSystemResto<OcpType>> &ls,
                                            VecRealView &x)
{
    const ProblemInfo<OcpType> &info = ls.info_;
    reduce(ls);
    VecRealView rhs_f_s = ls.rhs_f_s_.block(info.number_of_slack_variables, 0);
    VecRealView rhs_cl = ls.rhs_cl_.block(info.number_of_slack_variables, 0);
    VecRealView rhs_cu = ls.rhs_cu_.block(info.number_of_slack_variables, 0);

    LinearSystem<PdSystemType<OcpType>> ls_orig(
        ls.info_, ls.jac_, ls.hess_,
        ls.D_x_.block(info.number_of_primal_variables + info.number_of_slack_variables, 0), false,
        D_e_orig_, ls.Sl_i_.block(info.number_of_slack_variables, 0),
        ls.Su_i_.block(info.number_of_slack_variables, 0),
        ls.Zl_i_.block(info.number_of_slack_variables, 0),
        ls.Zu_i_.block(info.number_of_slack_variables, 0), ls.rhs_f_x_, rhs_f_s, rhs_g_orig_,
        rhs_cl, rhs_cu);
    // solve ls orig and save in x_orig_
    orig_solver_->solve_rhs(ls_orig, x_orig_);
    dereduce(ls, x);
}