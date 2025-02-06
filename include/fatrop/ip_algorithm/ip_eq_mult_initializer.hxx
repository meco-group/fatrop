//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_eq_mult_initializer_hxx__
#define __fatrop_ip_algorithm_ip_eq_mult_initializer_hxx__
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
          rhs_s_(ipdata->current_iterate().nlp()->nlp_dims().number_of_ineq_constraints),
          rhs_g_(ipdata->current_iterate().nlp()->nlp_dims().number_of_eq_constraints),
          rhs_cl_(ipdata->current_iterate().nlp()->nlp_dims().number_of_ineq_constraints),
          rhs_cu_(ipdata->current_iterate().nlp()->nlp_dims().number_of_ineq_constraints),
          Dx_(ipdata->current_iterate().nlp()->nlp_dims().number_of_variables +
              ipdata->current_iterate().nlp()->nlp_dims().number_of_ineq_constraints),
          Ds_(ipdata->current_iterate().nlp()->nlp_dims().number_of_ineq_constraints),
          Deq_(ipdata->current_iterate().nlp()->nlp_dims().number_of_eq_constraints),
          dummy_s_(ipdata->current_iterate().nlp()->nlp_dims().number_of_ineq_constraints),
          dummy_z_(ipdata->current_iterate().nlp()->nlp_dims().number_of_ineq_constraints)
    {
    }

    template <typename ProblemType> void IpEqMultInitializer<ProblemType>::initialize_eq_mult()
    {
        IpIterateType &curr_it = ipdata_->current_iterate();
        curr_it.set_dual_eq(VecRealScalar(curr_it.dual_eq().m(), 0.));
        rhs_x_.block(rhs_x_.m(), 0) = curr_it.dual_infeasibility_x();
        rhs_s_.block(rhs_s_.m(), 0) = curr_it.dual_infeasibility_s();
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
            curr_it.info(), curr_it.jacobian(), curr_it.zero_hessian(), Dx_, true, Deq_, dummy_s_,
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
        if (solved && norm_inf(curr_it.delta_dual_eq()) < lam_max_)
            curr_it.set_dual_eq(rhs_g_);
        else
        {
            curr_it.set_dual_eq(VecRealScalar(curr_it.dual_eq().m(), 0.));
        }
    }

    template <typename ProblemType> void IpEqMultInitializer<ProblemType>::reset()
    {
        // Empty implementation
    }
} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_eq_mult_initializer_hxx__
