//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#include "fatrop/ocp/pd_solver_orig.hpp"
#include "fatrop/linear_algebra/linear_solver.hxx"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/ocp/aug_system_solver.hpp"
#include "fatrop/ocp/problem_info.hpp"
using namespace fatrop;

// instantiate the template class
template class LinearSolver<PdSolverOrig<OcpType>, PdSystemType<OcpType>>;

PdSolverOrig<OcpType>::PdSolverOrig(const ProblemInfo<OcpType> &info,
                                    const std::shared_ptr<AugSystemSolver<OcpType>> &aug_system_solver)
    : LinearSolver<PdSolverOrig<OcpType>, PdSystemType<OcpType>>(
          LinearSystem<PdSystemType<OcpType>>::m(info)),
      sigma_inverse_(info.number_of_slack_variables), ss_(info.number_of_slack_variables),
      g_ii_(info.number_of_slack_variables), D_ii_(info.number_of_slack_variables),
      gg_(info.number_of_eq_constraints), x_aug_(info.number_of_primal_variables),
      mult_aug_(info.number_of_eq_constraints), aug_system_solver_(aug_system_solver)
{
}
void PdSolverOrig<OcpType>::reduce(LinearSystem<PdSystemType<OcpType>> &ls)
{
    VecRealView gi =
        ls.rhs_g_.block(ls.info_.number_of_slack_variables, ls.info_.offset_g_eq_slack);
    sigma_inverse_ =
        1. / (ls.D_x_.block(ls.info_.number_of_slack_variables, ls.info_.offset_slack) +
              1. / ls.Sl_i_ * ls.Zl_i_ + 1. / ls.Su_i_ * ls.Zu_i_);
    ss_ = ls.rhs_f_s_ + 1. / ls.Sl_i_ * ls.rhs_cl_ - 1. / ls.Su_i_ * ls.rhs_cu_;
    D_ii_ = sigma_inverse_ + ls.D_e_.block(ls.info_.number_of_g_eq_slack, ls.info_.offset_g_eq_slack);
    g_ii_ = gi + sigma_inverse_ * ss_;

    gg_.block(ls.info_.number_of_g_eq_path, ls.info_.offset_g_eq_path) =
        ls.rhs_g_.block(ls.info_.number_of_g_eq_path, ls.info_.offset_g_eq_path);
    gg_.block(ls.info_.number_of_g_eq_dyn, ls.info_.offset_g_eq_dyn) =
        ls.rhs_g_.block(ls.info_.number_of_g_eq_dyn, ls.info_.offset_g_eq_dyn);
    gg_.block(ls.info_.number_of_g_eq_slack, ls.info_.offset_g_eq_slack) = g_ii_;
}
void PdSolverOrig<OcpType>::dereduce(LinearSystem<PdSystemType<OcpType>> &ls, VecRealView &x)
{
    // set x
    x.block(ls.info_.number_of_primal_variables, ls.info_.pd_orig_offset_primal) = x_aug_;
    // set mult
    x.block(ls.info_.number_of_eq_constraints, ls.info_.pd_orig_offset_mult) = mult_aug_;
    // set s
    VecRealView mult_i =
        mult_aug_.block(ls.info_.number_of_slack_variables, ls.info_.offset_g_eq_slack);
    x.block(ls.info_.number_of_slack_variables, ls.info_.pd_orig_offset_slack) =
        sigma_inverse_ * (mult_i - ss_);
    // set zl and zu
    x.block(ls.info_.number_of_slack_variables, ls.info_.pd_orig_offset_zl) =
        1. / ls.Sl_i_ *
        (-ls.rhs_cl_ -
         ls.Zl_i_ * x.block(ls.info_.number_of_slack_variables, ls.info_.pd_orig_offset_slack));
    x.block(ls.info_.number_of_slack_variables, ls.info_.pd_orig_offset_zu) =
        1. / ls.Su_i_ *
        (-ls.rhs_cu_ +
         ls.Zu_i_ * x.block(ls.info_.number_of_slack_variables, ls.info_.pd_orig_offset_slack));
}
LinsolReturnFlag PdSolverOrig<OcpType>::solve_once_impl(LinearSystem<PdSystemType<OcpType>> &ls,
                                                        VecRealView &x)
{
    //    [ H + D_x    0        A_e^T  A_d^T  A_i^T    0     0  ] [ x   ] = [ -f_x ]
    //    [    0     D_x          0      0    -I      -I     I  ] [ s   ] = [ -f_s ]
    //    [ A_e       0        -D_e      0     0       0     0  ] [ λ_e ] = [ -g_e ]
    //    [ A_d       0          0       0     0       0     0  ] [ λ_d ] = [ -g_d ]
    //    [ A_i      -I          0       0   -D_i      0     0  ] [ λ_i ] = [ -g_i ]
    //    [   0     Zl_i         0       0     0     Sl_i    0  ] [ zl  ] = [ -cl  ]
    //    [   0    -Zu_i         0       0     0      0    Su_i ] [ zu  ] = [ -cu  ]

    // Step 1: Eliminate \(zl\) and \(zu\)
    //
    // From the last equation:
    //     \( zl = Sl_i^{-1} (-cl - Zl_i s) \)
    //     \( zu = Su_i^{-1} (-cu + Zu_i s) \)

    //    [ H + D_x    0        A_e^T  A_d^T  A_i^T ] [ x   ] = [ -f_x ]
    //    [    0     \Sigma      0       0    -I    ] [ s   ] = [ -ss  ]
    //    [ A_e       0        -D_e      0     0    ] [ λ_e ] = [ -g_e ]
    //    [ A_d       0          0       0     0    ] [ λ_d ] = [ -g_d ]
    //    [ A_i      -I          0       0    -D_i  ] [ λ_i ] = [ -g_i ]
    //    with \Sigma = D_x + Sl{-1} Zl + Su{-1} Zu
    //          ss = f_s  + Sl_i^{-1} cl - Su_i^{-1} cu
    // Step 2: Eliminate
    //    \(s\) = \Sigma{-1} (λ_i - ss)

    //    [ H + D_x    A_e^T  A_d^T  A_i^T  ] [ x   ] = [ -f_x ]
    //    [ A_e       -D_e      0     0     ] [ λ_e ] = [ -g_e ]
    //    [ A_d         0       0     0     ] [ λ_d ] = [ -g_d ]
    //    [ A_i         0       0     -D_ii ] [ λ_i ] = [ -g_ii]
    //    with D_ii = \Sigma^{-1} + D_i
    //         g_ii =  g_i + \Sigma{-1} ss
    //  This system is in the Augmented system form and can be solved by AugSystemSolver<OcpType>
    // call aug_system_solver to solve the system
    reduce(ls);
    LinsolReturnFlag ret;
    if (!ls.inertia_e_)
        ret = aug_system_solver_->solve(ls.info_, ls.jac_, ls.hess_, ls.D_x_, D_ii_, ls.rhs_f_x_,
                                        gg_, x_aug_, mult_aug_);
    else
        ret = aug_system_solver_->solve(ls.info_, ls.jac_, ls.hess_, ls.D_x_, ls.D_e_, D_ii_,
                                        ls.rhs_f_x_, gg_, x_aug_, mult_aug_);
    dereduce(ls, x);
    return ret;
}
void PdSolverOrig<OcpType>::solve_rhs_impl(LinearSystem<PdSystemType<OcpType>> &ls, VecRealView &x)
{
    reduce(ls);
    if (!ls.inertia_e_)
        aug_system_solver_->solve_rhs(ls.info_, ls.jac_, ls.hess_, D_ii_, ls.rhs_f_x_, gg_, x_aug_,
                                      mult_aug_);
    else
        aug_system_solver_->solve_rhs(ls.info_, ls.jac_, ls.hess_, ls.D_e_, D_ii_, ls.rhs_f_x_, gg_,
                                      x_aug_, mult_aug_);
    dereduce(ls, x);
}