//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_initializer_resto_hxx__
#define __fatrop_ip_algorithm_ip_initializer_resto_hxx__
#include "fatrop/common/options.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_eq_mult_initializer.hpp"
#include "fatrop/ip_algorithm/ip_initializer_resto.hpp"
#include "fatrop/nlp/nlp.hpp"

namespace fatrop
{
    template <typename ProblemType> void IpInitializerResto<ProblemType>::initialize()
    {
        IpDataType &data_orig = *data_orig_;
        IpDataType &data_resto = *data_resto_;
        const InfoType &info = data_orig.info();
        IpIterateType &curr_it_orig = data_orig.current_iterate();
        IpIterateType &curr_it_resto = data_resto.current_iterate();
        Index n_s = info.number_of_slack_variables;
        Index n_eq = info.number_of_eq_constraints;
        Index offset_n = info.offset_n;
        Index offset_p = info.offset_p;
        // compute barrier parameter
        Scalar resto_mu = std::max(curr_it_orig.mu(), norm_inf(curr_it_orig.constr_viol()));
        data_resto.current_iterate().set_mu(resto_mu);
        data_resto.trial_iterate().set_mu(resto_mu);
        // copy x from the original data to the resto data
        curr_it_resto.set_primal_x(curr_it_orig.primal_x());
        // copy s from the original data to the resto data
        // curr_it_resto.set_primal_s(curr_it_orig.primal_s().block(n_s, 0));
        curr_it_resto.primal_s().block(n_s, 0) = curr_it_orig.primal_s().block(n_s, 0);
        VecRealView constr = curr_it_orig.constr_viol();
        VecRealView p = curr_it_resto.primal_s().block(n_eq, offset_p);
        VecRealView n = curr_it_resto.primal_s().block(n_eq, offset_n);
        // xx =  (mu - rho* constr) / (2*rho)
        // p = xx * sqrt(xx**2 + mu*constr / (2rho))
        // n = constr + n
        p = (VecRealScalar(n_eq, resto_mu) - rho_ * constr) / VecRealScalar(n_eq, 2 * rho_);
        p = p + sqrt(p * p + resto_mu * constr / VecRealScalar(n_eq, 2 * rho_));
        n = constr + p;
        // intialize the bound multipliers
        VecRealView zl_orig = curr_it_orig.dual_bounds_l();
        VecRealView zu_orig = curr_it_orig.dual_bounds_u();
        VecRealView zl_resto = curr_it_resto.dual_bounds_l();
        VecRealView zu_resto = curr_it_resto.dual_bounds_u();

        zl_resto.block(n_s, 0) = min(zl_orig, VecRealScalar(n_s, rho_));
        zu_resto.block(n_s, 0) = min(zu_orig, VecRealScalar(n_s, rho_));
        zl_resto.block(n_eq, offset_p) = resto_mu * 1. / p;
        zl_resto.block(n_eq, offset_n) = resto_mu * 1. / n;
        // also set the upper bounds to 0. for convention
        zu_resto.block(n_eq, offset_p) = VecRealScalar(n_eq, 0.);
        zu_resto.block(n_eq, offset_n) = VecRealScalar(n_eq, 0.);
        // initialize the dual multipliers
        eq_mult_initializer_resto_->initialize_eq_mult();
    }

    template <typename ProblemType> void IpInitializerResto<ProblemType>::reset() {}

    template <typename ProblemType>
    void IpInitializerResto<ProblemType>::register_options(OptionRegistry &registry)
    {
        /**
         * todo register options
         */
    }

} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_initializer__resto_hxx__
