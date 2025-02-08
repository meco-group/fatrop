
//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_algorithm_hxx__
#define __fatrop_ip_algorithm_ip_algorithm_hxx__
#include "fatrop/ip_algorithm/ip_algorithm.hpp"
#include "fatrop/ip_algorithm/ip_convergence_check.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_eq_mult_initializer.hpp"
#include "fatrop/ip_algorithm/ip_initializer.hpp"
#include "fatrop/ip_algorithm/ip_iteration_output.hpp"
#include "fatrop/ip_algorithm/ip_linesearch.hpp"
#include "fatrop/ip_algorithm/ip_mu_update.hpp"
#include "fatrop/ip_algorithm/ip_search_dir.hpp"
#include "fatrop/ocp/type.hpp"
namespace fatrop
{
    template <typename ProblemType>
    IpAlgorithm<ProblemType>::IpAlgorithm(const IpSearchDirSp &search_dir,
                                          const IpLineSearchSp &linesearch,
                                          const IpInitializerSp &initializer,
                                          const IpMuUpdateSp &mu_update,
                                          const IpEqMultInitializerSp &eq_mult_initializer,
                                          const IpConvergenceCheckSp &convergence_check,
                                          const IpIterationOutputSp &iteration_output,
                                          const IpDataSp &ip_data)
        : search_dir_(search_dir), linesearch_(linesearch), initializer_(initializer),
          mu_update_(mu_update), eq_mult_initializer_(eq_mult_initializer),
          convergence_check_(convergence_check), iteration_output_(iteration_output),
          ip_data_(ip_data)
    {
    }

    template <typename ProblemType> void IpAlgorithm<ProblemType>::reset()
    {
        // todo who resets the ipdata?
        search_dir_->reset();
        linesearch_->reset();
        initializer_->reset();
        mu_update_->reset();
        eq_mult_initializer_->reset();
        convergence_check_->reset();
        // Note: IpIterationOutput might not need a reset method
    }

    template <typename ProblemType>
    IpSolverReturnFlag IpAlgorithm<ProblemType>::optimize(const bool is_resto)
    {
        reset();
        initializer_->initialize();
        IpSolverReturnFlag retval = IpSolverReturnFlag::Unknown;
        IpConvergenceStatus conv_status = convergence_check_->check_converged();

        iteration_output_->print_header();

        while (conv_status == IpConvergenceStatus::Continue)
        {
            iteration_output_->output_current_iteration();
            mu_update_->update_barrier_parameter();
            search_dir_->compute_search_dir();
            linesearch_->find_acceptable_trial_point();
            linesearch_->accept_trial_iterate();
            conv_status = convergence_check_->check_converged();
            ip_data_->set_iteration_number(ip_data_->iteration_number() + 1);
        }
        iteration_output_->output_current_iteration();
        if (conv_status == IpConvergenceStatus::Converged)
            retval = IpSolverReturnFlag::Success;
        return retval;
    }
    template <typename ProblemType>
    const ProblemInfo<ProblemType> &IpAlgorithm<ProblemType>::info() const
    {
        return ip_data_->current_iterate().info();
    }

    template <typename ProblemType>
    const VecRealView &IpAlgorithm<ProblemType>::solution_primal() const
    {
        return ip_data_->current_iterate().primal_x();
    }
    template <typename ProblemType>
    const VecRealView &IpAlgorithm<ProblemType>::solution_dual() const
    {
        return ip_data_->current_iterate().dual_eq();
    }


} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_data_hxx__
