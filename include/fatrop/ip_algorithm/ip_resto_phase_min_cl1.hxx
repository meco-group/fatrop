//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_resto_phase_min_cl1_hxx__
#define __fatrop_ip_resto_phase_min_cl1_hxx__
#include "fatrop/common/options.hpp"
#include "fatrop/ip_algorithm/ip_algorithm.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_iterate.hpp"
#include "fatrop/ip_algorithm/ip_resto_phase_min_cl1.hpp"
#include "fatrop/ip_algorithm/ip_eq_mult_initializer.hpp"

namespace fatrop
{
    template <typename ProblemType>
    IpRestoPhaseMinCl1<ProblemType>::IpRestoPhaseMinCl1(const IpAlgorithmSp &resto_ip_algorithm,
                                                        const IpInitializerBaseSp &ip_initializer,
                                                        const IpDataSp &ip_data_orig,
                                                        const IpDataSp &ip_data_resto,
                                                        const IpNlpRestoSp &ip_nlp_resto)
        : resto_ip_algorithm_(resto_ip_algorithm), ip_initializer_(ip_initializer),
          ip_data_orig_(ip_data_orig), ip_data_resto_(ip_data_resto), ip_nlp_resto_(ip_nlp_resto)
    {
    }

    template <typename ProblemType> bool IpRestoPhaseMinCl1<ProblemType>::perform_restoration()
    {
        count_resto_++;
        // initialize the ip_data_resto_ to have the data of the original problem
        ip_data_resto_->set_iteration_number(ip_data_orig_->iteration_number() + 1);
        IpIterateType &curr_it_resto = ip_data_resto_->current_iterate();
        IpIterateType &curr_it_orig = ip_data_orig_->current_iterate();
        // todo make sure that step info is set just before restoration phase is called in line
        // search
        curr_it_resto.step_info() = curr_it_orig.step_info();
        IpSolverReturnFlag resto_status = resto_ip_algorithm_->optimize(true);

        int retval = -1;

        IpIterateType &trial_it_resto = ip_data_resto_->trial_iterate();
        IpIterateType &trial_it_orig = ip_data_orig_->trial_iterate();
        const ProblemInfoType &info = curr_it_orig.info();

        if (resto_status != IpSolverReturnFlag::Success)
        {
            // set the current iterate of the original ip data such that it can be returned to the
            curr_it_orig.set_primal_x(curr_it_resto.primal_x());
            curr_it_orig.set_primal_s(
                curr_it_resto.primal_s().block(info.number_of_slack_varianbles, 0));
            // what to do with the duals?
            // Ipopt returns them but they dont really make sense for the original problem.
            fatrop_assert(false && "Restoration phase failed");
        }
        if (resto_status == IpSolverReturnFlag::Success)
        {
            retval = 0;
        }
        if (retval == 0)
        {
            // compute delta s from curr_s of the original solver to curr_s of the resto solver,
            // thsi will be used later to compute delta z
            trial_it_orig.set_delta_primal_s(
                curr_it_resto.primal_s().block(info.number_of_slack_varianbles, 0) -
                curr_it_orig.primal_s().block(info.number_of_slack_varianbles, 0));
            // now we compute delta z as if one step was taken.
            const VecRealView Sl = curr_it_orig.delta_lower();
            const VecRealView Su = curr_it_orig.delta_upper();
            const VecRealView Zl = curr_it_orig.delta_dual_bounds_l();
            const VecRealView Zu = curr_it_orig.delta_dual_bounds_u();
            const VecRealView compl_l = curr_it_orig.relaxed_complementarity_l();
            const VecRealView compl_u = curr_it_orig.relaxed_complementarity_u();
            const VecRealView delta_s = trial_it_orig.delta_primal_s();
            trial_it_orig.set_delta_dual_bounds_l(1. / Sl * (-compl_l - Zl * delta_s));
            trial_it_orig.set_delta_dual_bounds_u(1. / Su * (-compl_u + Zu * delta_s));
            // now compute the max step size
            Scalar alpha_dual = curr_it_orig.maximum_step_size_dual(
                curr_it_orig.tau(), trial_it_orig.dual_bounds_l(), trial_it_orig.dual_bounds_u());
            // set lower and upper bounds for the step size
            trial_it_orig.set_dual_bounds_l(curr_it_orig.dual_bounds_l() +
                                            alpha_dual * trial_it_orig.delta_dual_bounds_l());
            trial_it_orig.set_dual_bounds_u(curr_it_orig.dual_bounds_u() +
                                            alpha_dual * trial_it_orig.delta_dual_bounds_u());
            Scalar bound_mult_max = std::max(norm_inf(trial_it_orig.dual_bounds_l()),
                                             norm_inf(trial_it_orig.dual_bounds_u()));
            if (bound_mult_max > bound_mult_reset_treshold_)
            {
                trial_it_orig.set_dual_bounds_l(
                    VecRealScalar(trial_it_orig.dual_bounds_l().m(), 1.));
                trial_it_orig.set_dual_bounds_u(
                    VecRealScalar(trial_it_orig.dual_bounds_u().m(), 1.));
            }

            // copy the result to the trial iterate of the original problem
            trial_it_orig.set_primal_x(curr_it_resto.primal_x());
            trial_it_orig.set_primal_s(
                curr_it_resto.primal_s().block(info.number_of_slack_varianbles, 0));
            // run the eq mult initializer
            ip_initializer_->initialize_eq_mult(true);
            ip_data_orig_->set_iteration_number(ip_data_orig_->iteration_number() - 1);
        }

        return (retval == 0);
    }

    template <typename ProblemType> void IpRestoPhaseMinCl1<ProblemType>::reset()
    {
        count_resto_ = 0;
    }

    template <typename ProblemType>
    void IpRestoPhaseMinCl1<ProblemType>::register_options(OptionRegistry &registry)
    {
        registry.register_option("bound_mult_reset_treshold",
                                 &IpRestoPhaseMinCl1::set_bound_mult_reset_treshold, this);
        registry.register_option("constr_mult_reset_treshold",
                                 &IpRestoPhaseMinCl1::set_constr_mult_reset_treshold, this);
        registry.register_option("resto_failure_feasibility_treshold",
                                 &IpRestoPhaseMinCl1::set_resto_failure_feasibility_treshold, this);
        registry.register_option("constr_viol_tol", &IpRestoPhaseMinCl1::set_constr_viol_tol, this);
    }

} // namespace fatrop

#endif // __fatrop_ip_resto_phase_min_cl1_hxx__
