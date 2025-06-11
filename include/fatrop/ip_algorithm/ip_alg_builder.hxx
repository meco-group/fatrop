#ifndef __fatrop_ip_algorithm_ip_alg_builder_hxx__
#define __fatrop_ip_algorithm_ip_alg_builder_hxx__

#include "fatrop/common/options.hpp"
#include "fatrop/ip_algorithm/ip_alg_builder.hpp"
#include "fatrop/ip_algorithm/ip_algorithm.hpp"
#include "fatrop/ip_algorithm/ip_convergence_check.hpp"
#include "fatrop/ip_algorithm/ip_convergence_check_resto.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_eq_mult_initializer.hpp"
#include "fatrop/ip_algorithm/ip_initializer.hpp"
#include "fatrop/ip_algorithm/ip_initializer_resto.hpp"
#include "fatrop/ip_algorithm/ip_iteration_output.hpp"
#include "fatrop/ip_algorithm/ip_iteration_output_resto.hpp"
#include "fatrop/ip_algorithm/ip_linesearch.hpp"
#include "fatrop/ip_algorithm/ip_mu_update.hpp"
#include "fatrop/ip_algorithm/ip_nlp_orig.hpp"
#include "fatrop/ip_algorithm/ip_nlp_resto.hpp"
#include "fatrop/ip_algorithm/ip_resto_phase_min_cl1.hpp"
#include "fatrop/ip_algorithm/ip_search_dir.hpp"
#include "fatrop/ip_algorithm/ip_timings.hpp"
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
        /**
         * todo problem info also gets created by the constructor o fthe IpData
         */
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
        search_dir_ = std::make_shared<IpSearchDirImpl<PdSolverOrig<ProblemType>, ProblemType>>(
            ipdata_, pd_solver_);
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
        if (!resto_phase_)
            create_restoration_phase();
        linesearch_ = std::make_shared<IpLinesearch<PdSolverOrig<ProblemType>, ProblemType>>(
            ipdata_, pd_solver_, resto_phase_);
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
    IpAlgBuilder<ProblemType> &IpAlgBuilder<ProblemType>::create_restoration_phase()
    {
        /**
         *
         * create the algorithm for the restoration phase
         *
         */
        if (!problem_info_)
            create_problem_info();
        if (!pd_solver_)
            create_pdsolver();
        if (!ipdata_)
            create_ipdata();
        if (!iteration_output_)
            create_iteration_output();
        // create the resto nlp
        std::shared_ptr<IpNlpResto<ProblemType>> ip_nlp_resto =
            std::make_shared<IpNlpResto<ProblemType>>(nlp_orig_);
        if (options_registry_)
            options_registry_->register_options(*ip_nlp_resto);
        // create the resto ipdata
        std::shared_ptr<IpData<ProblemType>> ip_data_resto =
            std::make_shared<IpData<ProblemType>>(ip_nlp_resto);
        if (options_registry_)
            options_registry_->register_options(*ip_data_resto);
        // create the PdSolverResto
        std::shared_ptr<PdSolverResto<ProblemType>> pd_solver_resto =
            std::make_shared<PdSolverResto<ProblemType>>(*problem_info_, pd_solver_);
        // create the resto search dir
        typedef IpSearchDirImpl<PdSolverResto<ProblemType>, ProblemType> RestoSDType;
        std::shared_ptr<RestoSDType> search_dir_resto =
            std::make_shared<RestoSDType>(ip_data_resto, pd_solver_resto);
        if (options_registry_)
            options_registry_->register_options(*search_dir_resto);
        // create the resto linesearch
        typedef IpLinesearch<PdSolverResto<ProblemType>, ProblemType> RestoLSType;
        std::shared_ptr<RestoLSType> linesearch_resto =
            std::make_shared<RestoLSType>(ip_data_resto, pd_solver_resto, nullptr);
        if (options_registry_)
            options_registry_->register_options(*linesearch_resto);
        // create the convergence check for the resto phase
        convergence_check_resto_ =
            std::make_shared<IpConvergenceCheckResto<ProblemType>>(ipdata_, ip_data_resto);
        if (options_registry_)
            options_registry_->register_options(*convergence_check_resto_);
        // create the mu update for the resto phase
        std::shared_ptr<IpMonotoneMuUpdate<ProblemType>> mu_update_resto =
            std::make_shared<IpMonotoneMuUpdate<ProblemType>>(ip_data_resto, linesearch_resto);
        if (options_registry_)
            options_registry_->register_options(*mu_update_resto);
        // create eq mult initializer
        std::shared_ptr<IpEqMultInitializer<ProblemType>> eq_mult_initializer_resto =
            std::make_shared<IpEqMultInitializer<ProblemType>>(ip_data_resto, pd_solver_);
        if (options_registry_)
            options_registry_->register_options(*eq_mult_initializer_resto);
        // create the resto initializer
        std::shared_ptr<IpInitializerResto<ProblemType>> initializer_resto =
            std::make_shared<IpInitializerResto<ProblemType>>(ipdata_, ip_data_resto);
        if (options_registry_)
            options_registry_->register_options(*initializer_resto);
        // create iteration output for the resto phase
        std::shared_ptr<IpIterationOutputResto<ProblemType>> iteration_output_resto =
            std::make_shared<IpIterationOutputResto<ProblemType>>(ipdata_, ip_data_resto);
        if (options_registry_)
            options_registry_->register_options(*iteration_output_resto);
        // create the resto algorithn
        std::shared_ptr<IpAlgorithm<ProblemType>> resto_alg =
            std::make_shared<IpAlgorithm<ProblemType>>(
                search_dir_resto, linesearch_resto, initializer_resto, mu_update_resto,
                eq_mult_initializer_resto, convergence_check_resto_, iteration_output_resto,
                ip_data_resto);
        resto_phase_ = std::make_shared<IpRestoPhaseMinCl1<ProblemType>>(
            resto_alg, eq_mult_initializer_resto, ipdata_, ip_data_resto, ip_nlp_resto);
        if (options_registry_)
            options_registry_->register_options(*resto_phase_);
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
        if (!resto_phase_)
            create_restoration_phase();
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
        // set the timings satistics
        nlp_orig_->set_timing_statistics(&ipdata_->timing_statistics());
        if (resto_phase_)
            convergence_check_resto_->set_line_search_orig(linesearch_);

        return std::make_shared<IpAlgorithm<ProblemType>>(
            search_dir_, linesearch_, initializer_, mu_update_, eq_mult_initializer_,
            convergence_check_, iteration_output_, ipdata_);
    }
} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_alg_builder_hxx__
