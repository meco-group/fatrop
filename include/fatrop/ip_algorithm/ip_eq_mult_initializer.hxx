//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_eq_mult_initializer_hxx__
#define __fatrop_ip_algorithm_ip_eq_mult_initializer_hxx__
#include "fatrop/common/options.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_eq_mult_initializer.hpp"
#include "fatrop/ip_algorithm/pd_solver_orig.hpp"
#include "fatrop/ip_algorithm/pd_system_orig.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"

namespace fatrop
{
    template <typename ProblemType>
    IpEqMultInitializer<ProblemType>::IpEqMultInitializer(const IpDataSp &ipdata,
                                                          const PdSolverSp &linear_solver)
        : ipdata_(ipdata), linear_solver_(linear_solver),
          rhs_x_(ipdata->current_iterate().nlp()->nlp_dims().number_of_variables),
          rhs_s_(ipdata->current_iterate().info().number_of_slack_variables),
          rhs_g_(ipdata->current_iterate().nlp()->nlp_dims().number_of_eq_constraints),
          rhs_cl_(ipdata->current_iterate().info().number_of_slack_variables),
          rhs_cu_(ipdata->current_iterate().info().number_of_slack_variables),
          Dx_(ipdata->current_iterate().nlp()->nlp_dims().number_of_variables +
              ipdata->current_iterate().info().number_of_slack_variables),
          Ds_(ipdata->current_iterate().info().number_of_slack_variables),
          Deq_(ipdata->current_iterate().nlp()->nlp_dims().number_of_eq_constraints),
          dummy_s_(ipdata->current_iterate().info().number_of_slack_variables),
          dummy_z_(ipdata->current_iterate().info().number_of_slack_variables)
    {
    }

    template <typename ProblemType>
    void IpEqMultInitializer<ProblemType>::initialize_eq_mult(bool trial_it /* = false */)
    {
        IpIterateType &iterate = trial_it ? ipdata_->trial_iterate() : ipdata_->current_iterate();
        iterate.set_dual_eq(VecRealScalar(iterate.dual_eq().m(), 0.));
        rhs_x_.block(rhs_x_.m(), 0) = iterate.dual_infeasibility_x();
        rhs_s_.block(rhs_s_.m(), 0) = iterate.dual_infeasibility_s().block(
            ipdata_->current_iterate().info().number_of_slack_variables, 0);
        rhs_g_.block(rhs_g_.m(), 0) = 0;
        rhs_cl_.block(rhs_cl_.m(), 0) = 0;
        rhs_cu_.block(rhs_cu_.m(), 0) = 0;

        Dx_ = 1.;
        Deq_ = 0.;
        dummy_s_ = 1.;
        dummy_z_ = 0.;

        bool solved = false;
        bool first_try_delta_w = true;
        LinearSystem<PdSystemType<ProblemType>> ls(
            iterate.info(), iterate.jacobian(), iterate.zero_hessian(), Dx_, true, Deq_, dummy_s_,
            dummy_s_, dummy_z_, dummy_z_, rhs_x_, rhs_s_, rhs_g_, rhs_cl_, rhs_cu_);
        LinsolReturnFlag ret = linear_solver_->solve_in_place(ls);
        switch (ret)
        {
        case (LinsolReturnFlag::SUCCESS):
            solved = true;
            break;
        case (LinsolReturnFlag::ITREF_MAX_ITER):
            solved = true;
            break;
        case (LinsolReturnFlag::ITREF_INCREASE):
            solved = true;
            break;
        case (LinsolReturnFlag::INDEFINITE):
            solved = false;
            break;
        case (LinsolReturnFlag::NOFULL_RANK):
            solved = false;
            break;
        case (LinsolReturnFlag::UNKNOWN):
            fatrop_assert_msg(false, "Unexpected return flag from linear solver");
            break;
        }
        if (solved && norm_inf(rhs_g_) < lam_max_)
            if (trial_it)
                iterate.set_dual_eq(rhs_g_);
            else
                iterate.set_dual_eq(rhs_g_);
        else
        {
            iterate.set_dual_eq(VecRealScalar(iterate.dual_eq().m(), 0.));
        }
    }

    template <typename ProblemType> void IpEqMultInitializer<ProblemType>::reset()
    {
        // Empty implementation
    }
    template <typename ProblemType>
    void IpEqMultInitializer<ProblemType>::register_options(OptionRegistry &registry)
    {
        registry.register_option("lam_max", &IpEqMultInitializer::set_lam_max, this);
    }

} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_eq_mult_initializer_hxx__
