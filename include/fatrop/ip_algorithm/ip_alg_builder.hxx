#ifndef __fatrop_ip_algorithm_ip_alg_builder_hxx__
#define __fatrop_ip_algorithm_ip_alg_builder_hxx__

#include "fatrop/common/options.hpp"
#include "fatrop/ip_algorithm/ip_alg_builder.hpp"
#include "fatrop/ip_algorithm/ip_algorithm.hpp"
#include "fatrop/ip_algorithm/ip_convergence_check.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_eq_mult_initializer.hpp"
#include "fatrop/ip_algorithm/ip_initializer.hpp"
#include "fatrop/ip_algorithm/ip_iteration_output.hpp"
#include "fatrop/ip_algorithm/ip_linesearch.hpp"
#include "fatrop/ip_algorithm/ip_mu_update.hpp"
#include "fatrop/ip_algorithm/ip_nlp_orig.hpp"
#include "fatrop/ip_algorithm/ip_search_dir.hpp"
#include "fatrop/linear_algebra/linear_solver.hpp"

namespace fatrop
{
    template <typename ProblemType>
    IpAlgBuilder<ProblemType>::IpAlgBuilder(const std::shared_ptr<Nlp<ProblemType>> &nlp)
    {
        nlp_orig_ = std::make_shared<IpNlpOrig<ProblemType>>(nlp);
    }

    template <typename ProblemType>
    IpAlgBuilder<ProblemType> &IpAlgBuilder<ProblemType>::create_iteration_output()
    {
        if (!ipdata_)
            create_ipdata();
        iteration_output_ = std::make_shared<IpIterationOutput<ProblemType>>(ipdata_);
        return *this;
    }

    template <typename ProblemType>
    IpAlgBuilder<ProblemType> &IpAlgBuilder<ProblemType>::create_ipdata()
    {
        ipdata_ = std::make_shared<IpData<ProblemType>>(nlp_orig_);
        if (options_registry_)
            ipdata_->register_options(*options_registry_);
        return *this;
    }

    template <typename ProblemType>
    IpAlgBuilder<ProblemType> &IpAlgBuilder<ProblemType>::create_problem_info()
    {
        const auto &ocp_dims = nlp_orig_->problem_dims();
        problem_info_ = std::make_shared<ProblemInfo<ProblemType>>(ocp_dims);
        return *this;
    }

    template <typename ProblemType>
    IpAlgBuilder<ProblemType> &IpAlgBuilder<ProblemType>::create_aug_system_solver()
    {
        if (!problem_info_)
            create_problem_info();
        aug_system_solver_ = std::make_shared<AugSystemSolver<ProblemType>>(*problem_info_);
        return *this;
    }

    template <typename ProblemType>
    IpAlgBuilder<ProblemType> &IpAlgBuilder<ProblemType>::create_pdsolver()
    {
        if (!problem_info_)
            create_problem_info();
        if (!aug_system_solver_)
            create_aug_system_solver();
        pd_solver_ =
            std::make_shared<PdSolverOrig<ProblemType>>(*problem_info_, aug_system_solver_);
        return *this;
    }

    template <typename ProblemType>
    IpAlgBuilder<ProblemType> &IpAlgBuilder<ProblemType>::create_search_dir()
    {
        if (!ipdata_)
            create_ipdata();
        if (!pd_solver_)
            create_pdsolver();
        search_dir_ = std::make_shared<IpSearchDirImpl<PdSolverOrig<ProblemType>, ProblemType>>(ipdata_, pd_solver_);
        if (options_registry_)
            options_registry_->register_options(*search_dir_);
        return *this;
    }

    template <typename ProblemType>
    IpAlgBuilder<ProblemType> &IpAlgBuilder<ProblemType>::create_linesearch()
    {
        if (!ipdata_)
            create_ipdata();
        if (!pd_solver_)
            create_pdsolver();
        linesearch_ = std::make_shared<IpLinesearch<PdSolverOrig<ProblemType>, ProblemType>>(ipdata_, pd_solver_);
        if (options_registry_)
            options_registry_->register_options(*linesearch_);
        return *this;
    }

    template <typename ProblemType>
    IpAlgBuilder<ProblemType> &IpAlgBuilder<ProblemType>::create_initializer()
    {
        if (!ipdata_)
            create_ipdata();
        if (!eq_mult_initializer_)
            create_eq_mult_initializer();
        initializer_ = std::make_shared<IpInitializer<ProblemType>>(ipdata_, eq_mult_initializer_);
        if (options_registry_)
            options_registry_->register_options(*initializer_);
        return *this;
    }

    template <typename ProblemType>
    IpAlgBuilder<ProblemType> &IpAlgBuilder<ProblemType>::create_mu_update()
    {
        if (!ipdata_)
            create_ipdata();
        if (!linesearch_)
            create_linesearch();
        mu_update_ = std::make_shared<IpMonotoneMuUpdate<ProblemType>>(ipdata_, linesearch_);
        if (options_registry_)
            options_registry_->register_options(*mu_update_);
        return *this;
    }

    template <typename ProblemType>
    IpAlgBuilder<ProblemType> &IpAlgBuilder<ProblemType>::create_eq_mult_initializer()
    {
        if (!ipdata_)
            create_ipdata();
        if (!pd_solver_)
            create_pdsolver();
        eq_mult_initializer_ =
            std::make_shared<IpEqMultInitializer<ProblemType>>(ipdata_, pd_solver_);
        if (options_registry_)
            options_registry_->register_options(*eq_mult_initializer_);
        return *this;
    }

    template <typename ProblemType>
    IpAlgBuilder<ProblemType> &IpAlgBuilder<ProblemType>::create_convergence_check()
    {
        if (!ipdata_)
            create_ipdata();
        convergence_check_ = std::make_shared<IpConvergenceCheck<ProblemType>>(ipdata_);
        if (options_registry_)
            options_registry_->register_options(*convergence_check_);
        return *this;
    }

    template <typename ProblemType>
    std::shared_ptr<IpAlgorithm<ProblemType>> IpAlgBuilder<ProblemType>::build()
    {
        // if (!ipdata_)
        //     create_ipdata();
        // if (!problem_info_)
        //     create_problem_info();
        // if (!aug_system_solver_)
        //     create_aug_system_solver();
        // if (!pd_solver_)
        //     create_pdsolver();
        if (!search_dir_)
            create_search_dir();
        if (!linesearch_)
            create_linesearch();
        if (!initializer_)
            create_initializer();
        if (!mu_update_)
            create_mu_update();
        if (!eq_mult_initializer_)
            create_eq_mult_initializer();
        if (!convergence_check_)
            create_convergence_check();
        if (!iteration_output_)
            create_iteration_output();

        return std::make_shared<IpAlgorithm<ProblemType>>(
            search_dir_, linesearch_, initializer_, mu_update_, eq_mult_initializer_,
            convergence_check_, iteration_output_, ipdata_);
    }
} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_alg_builder_hxx__
