//
// Copyright (c) 2026 Lander Vanroye, KU Leuven
//

#include "fatrop/graph/pd_solver_orig.hpp"

#include "fatrop/graph/aug_system_solver.hpp"
#include "fatrop/graph/problem_info.hpp"
#include "fatrop/linear_algebra/linear_solver.hxx"
#include "fatrop/linear_algebra/vector.hpp"

using namespace fatrop;

// Instantiate the CRTP base template for GraphProblem.
template class fatrop::LinearSolver<PdSolverOrig<GraphProblem>, PdSystemType<GraphProblem>>;

PdSolverOrig<GraphProblem>::PdSolverOrig(
    const ProblemInfo<GraphProblem> &info,
    const std::shared_ptr<AugSystemSolver<GraphProblem>> &aug_system_solver)
    : LinearSolver<PdSolverOrig<GraphProblem>, PdSystemType<GraphProblem>>(
          LinearSystem<PdSystemType<GraphProblem>>::m(info)),
      sigma_inverse_(info.number_of_slack_variables), ss_(info.number_of_slack_variables),
      g_ii_(info.number_of_slack_variables), D_ii_(info.number_of_slack_variables),
      gg_(info.number_of_eq_constraints), x_aug_(info.number_of_tangent_variables),
      mult_aug_(info.number_of_eq_constraints), aug_system_solver_(aug_system_solver)
{
}

void PdSolverOrig<GraphProblem>::reduce(LinearSystem<PdSystemType<GraphProblem>> &ls)
{
    // For graph problems offset_g_eq_slack == 0 and number_of_g_eq_slack ==
    // number_of_slack_variables, so the slack-eq block coincides with the
    // entire equality vector.
    VecRealView gi = ls.rhs_g_.block(ls.info_.number_of_slack_variables, ls.info_.offset_g_eq_slack);
    sigma_inverse_ =
        1. / (ls.D_x_.block(ls.info_.number_of_slack_variables, ls.info_.offset_slack) +
              1. / ls.Sl_i_ * ls.Zl_i_ + 1. / ls.Su_i_ * ls.Zu_i_);
    ss_ = ls.rhs_f_s_ + 1. / ls.Sl_i_ * ls.rhs_cl_ - 1. / ls.Su_i_ * ls.rhs_cu_;
    D_ii_ = sigma_inverse_ + ls.D_e_.block(ls.info_.number_of_g_eq_slack, ls.info_.offset_g_eq_slack);
    g_ii_ = gi + sigma_inverse_ * ss_;

    // No path-equality block: gg_ == g_ii_ (number_of_eq_constraints ==
    // number_of_g_eq_slack).
    gg_.block(ls.info_.number_of_g_eq_slack, ls.info_.offset_g_eq_slack) = g_ii_;
}

void PdSolverOrig<GraphProblem>::dereduce(LinearSystem<PdSystemType<GraphProblem>> &ls,
                                           VecRealView &x)
{
    x.block(ls.info_.number_of_tangent_variables, ls.info_.pd_orig_offset_primal) = x_aug_;
    x.block(ls.info_.number_of_eq_constraints, ls.info_.pd_orig_offset_mult) = mult_aug_;
    VecRealView mult_i =
        mult_aug_.block(ls.info_.number_of_slack_variables, ls.info_.offset_g_eq_slack);
    x.block(ls.info_.number_of_slack_variables, ls.info_.pd_orig_offset_slack) =
        sigma_inverse_ * (mult_i - ss_);
    x.block(ls.info_.number_of_slack_variables, ls.info_.pd_orig_offset_zl) =
        1. / ls.Sl_i_ *
        (-ls.rhs_cl_ -
         ls.Zl_i_ * x.block(ls.info_.number_of_slack_variables, ls.info_.pd_orig_offset_slack));
    x.block(ls.info_.number_of_slack_variables, ls.info_.pd_orig_offset_zu) =
        1. / ls.Su_i_ *
        (-ls.rhs_cu_ +
         ls.Zu_i_ * x.block(ls.info_.number_of_slack_variables, ls.info_.pd_orig_offset_slack));
}

LinsolReturnFlag PdSolverOrig<GraphProblem>::solve_once_impl(
    LinearSystem<PdSystemType<GraphProblem>> &ls, VecRealView &x)
{
    reduce(ls);
    LinsolReturnFlag ret = aug_system_solver_->solve(ls.info_, ls.jac_, ls.hess_, ls.D_x_, D_ii_,
                                                     ls.rhs_f_x_, gg_, x_aug_, mult_aug_);
    dereduce(ls, x);
    return ret;
}

void PdSolverOrig<GraphProblem>::solve_rhs_impl(LinearSystem<PdSystemType<GraphProblem>> &ls,
                                                 VecRealView &x)
{
    reduce(ls);
    aug_system_solver_->solve_rhs(ls.info_, ls.jac_, ls.hess_, D_ii_, ls.rhs_f_x_, gg_, x_aug_,
                                  mult_aug_);
    dereduce(ls, x);
}
