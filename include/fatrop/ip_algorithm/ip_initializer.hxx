//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_initializer_hxx__
#define __fatrop_ip_algorithm_ip_initializer_hxx__
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_eq_mult_initializer.hpp"
#include "fatrop/ip_algorithm/ip_initializer.hpp"
#include "fatrop/nlp/nlp.hpp"
namespace fatrop
{

    template <typename ProblemType>
    IpInitializer<ProblemType>::IpInitializer(const IpDataSp ipdata,
                                              const IpEqMultInitializerSp &eq_mult_initializer)
        : ipdata_(ipdata), eq_mult_initializer_(eq_mult_initializer)
    {
    }

    template <typename ProblemType> void IpInitializer<ProblemType>::initialize()
    {
        const Index m = ipdata_->current_iterate().primal_s().m();
        const std::vector<bool> lower_bounded = ipdata_->current_iterate().lower_bounded();
        const std::vector<bool> upper_bounded = ipdata_->current_iterate().upper_bounded();
        // set z to 1. if bounded and 0. otherwise
        ipdata_->current_iterate().set_dual_bounds_l(if_else(lower_bounded, VecRealScalar(m, 1.0), VecRealScalar(m, 0.0)));
        ipdata_->current_iterate().set_dual_bounds_u(if_else(upper_bounded, VecRealScalar(m, 1.0), VecRealScalar(m, 0.0)));
        initialize_slacks();
        eq_mult_initializer_->initialize_eq_mult();
    }
    template <typename ProblemType> void IpInitializer<ProblemType>::initialize_slacks()
    {
        fatrop_assert_msg(ipdata_->current_iterate().primal_s().is_zero(),
                          "Slack variables must be zero at initialization");
        const VecRealView viol_s = ipdata_->current_iterate().constr_viol_ineq();
        const VecRealView lower_bounds = ipdata_->current_iterate().lower_bounds();
        const VecRealView upper_bounds = ipdata_->current_iterate().upper_bounds();
        const std::vector<bool> lower_bounded = ipdata_->current_iterate().lower_bounded();
        const std::vector<bool> upper_bounded = ipdata_->current_iterate().upper_bounded();
        auto double_bounded = [&](Index i) { return lower_bounded[i] && upper_bounded[i]; };
        auto single_bounded = [&](Index i) { return lower_bounded[i] ^ upper_bounded[i]; };
        Index number_of_slacks = lower_bounds.m();
        auto pl = min(bound_push * max(VecRealScalar(number_of_slacks, 1.0), abs(lower_bounds)),
                      bound_frac * (upper_bounds - lower_bounds));
        auto pr = min(bound_push * max(VecRealScalar(number_of_slacks, 1.0), abs(upper_bounds)),
                      bound_frac * (upper_bounds - lower_bounds));
        //   project the slacks
        auto res_double_bounded = min(max(viol_s, lower_bounds + pl), upper_bounds - pr);
        auto res_single_lower_bounded =
            max(viol_s, lower_bounds + bound_push * max(VecRealScalar(number_of_slacks, 1.0),
                                                        abs(lower_bounds)));
        auto res_single_upper_bounded =
            min(viol_s, upper_bounds - bound_push * max(VecRealScalar(number_of_slacks, 1.0),
                                                        abs(upper_bounds)));
        auto res =
            if_else(double_bounded, res_double_bounded,
                    if_else(lower_bounded, res_single_lower_bounded, res_single_upper_bounded));
        ipdata_->current_iterate().set_primal_s(res);
    }

    template <typename ProblemType> void IpInitializer<ProblemType>::reset()
    {
        // Empty implementation
    }
} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_initializer_hxx__
