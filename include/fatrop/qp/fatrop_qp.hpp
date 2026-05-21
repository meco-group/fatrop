//
// Copyright (C) 2026 Lander Vanroye, KU Leuven
//
// Public API for the Mehrotra predictor-corrector QP solver. The usage
// mirrors the NLP solver (IpAlgBuilder<OcpType> / IpAlgorithm<OcpType>):
//
//     OptionRegistry options;
//     MehrotraQpBuilder<OcpType> builder(std::make_shared<NlpOcp>(ocp));
//     auto qpalg = builder.with_options_registry(&options).build();
//     IpSolverReturnFlag ret = qpalg->optimize();
//     auto data = builder.get_ipdata();
//     std::cout << data->timing_statistics();
//
// The exit-status enum, options registry, and IpData / IpAlgorithm types are
// shared with the NLP solver -- the QP module just plugs in a different
// search-direction / line-search / mu-update strategy.
//
#ifndef __fatrop_qp_fatrop_qp_hpp__
#define __fatrop_qp_fatrop_qp_hpp__

#include "fatrop/common/fwd.hpp"
#include "fatrop/context/context.hpp"
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/ip_algorithm/ip_algorithm.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/nlp/fwd.hpp"

#include <memory>

namespace fatrop
{
    // ip_algorithm/fwd.hpp does not forward-declare IpInitializer.
    template <typename ProblemType> class IpInitializer;

    /**
     * @brief Mehrotra predictor-corrector QP solver for OCP-structured QPs.
     *
     * Same role as @c IpAlgorithm for the nonlinear solver, but each
     * iteration uses one factorization + two back-substitutions (predictor
     * and corrector) and a fraction-to-boundary step instead of a filter
     * line search.
     *
     * @tparam ProblemType Problem tag (e.g. @c OcpType).
     */
    template <typename ProblemType> class MehrotraQpAlgorithm
    {
    public:
        typedef std::shared_ptr<Nlp<ProblemType>> NlpSp;
        typedef std::shared_ptr<IpData<ProblemType>> IpDataSp;
        typedef std::shared_ptr<IpNlpOrig<ProblemType>> IpNlpOrigSp;
        typedef std::shared_ptr<PdSolverOrig<ProblemType>> PdSolverSp;
        typedef std::shared_ptr<IpInitializer<ProblemType>> InitializerSp;
        typedef std::shared_ptr<IpEqMultInitializer<ProblemType>> EqMultInitSp;

        MehrotraQpAlgorithm(const IpDataSp &ipdata, const IpNlpOrigSp &nlp_orig,
                            const PdSolverSp &pd_solver, const InitializerSp &initializer,
                            const EqMultInitSp &eq_mult_initializer);

        /**
         * @brief Run the predictor-corrector loop until convergence.
         */
        IpSolverReturnFlag optimize();

        /**
         * @brief Reset internal state (called automatically by @c optimize).
         */
        void reset();

        // -------- Options ---------------------------------------------------
        void set_tau_min(const Scalar &v) { tau_min_ = v; }
        void set_max_iter(const Index &v) { max_iter_ = v; }
        void set_tolerance(const Scalar &v) { tol_ = v; }
        void set_constr_viol_tol(const Scalar &v) { constr_viol_tol_ = v; }
        void set_mu_min(const Scalar &v) { mu_min_ = v; }
        void set_sigma_max(const Scalar &v) { sigma_max_ = v; }
        void set_verbose(bool v) { verbose_ = v; }

        // -------- Accessors -------------------------------------------------
        const ProblemInfo<ProblemType> &info() const;
        const VecRealView &solution_primal() const;
        const VecRealView &solution_dual() const;
        Index iteration_count() const { return iteration_; }
        Scalar final_mu() const { return mu_; }

        void register_options(OptionRegistry &registry);

    private:
        LinsolReturnFlag compute_predictor_corrector_step(Scalar &alpha_pr_out,
                                                          Scalar &alpha_du_out,
                                                          Scalar &sigma_out);
        void apply_step(Scalar alpha_pr, Scalar alpha_du);
        void recompute_mu();
        Scalar compute_mu_aff(Scalar alpha_pr, Scalar alpha_du, const VecRealView &ds_aff,
                              const VecRealView &dzl_aff, const VecRealView &dzu_aff) const;
        void print_header() const;
        void print_iteration(Scalar inf_pr, Scalar inf_du, Scalar inf_compl, Scalar mu_aff,
                             Scalar sigma, Scalar alpha_pr, Scalar alpha_du) const;

        IpDataSp ipdata_;
        IpNlpOrigSp nlp_orig_;
        PdSolverSp pd_solver_;
        InitializerSp initializer_;
        EqMultInitSp eq_mult_initializer_;

        VecRealAllocated rhs_x_;
        VecRealAllocated rhs_s_;
        VecRealAllocated rhs_g_;
        VecRealAllocated rhs_cl_;
        VecRealAllocated rhs_cu_;
        VecRealAllocated ds_aff_;
        VecRealAllocated dzl_aff_;
        VecRealAllocated dzu_aff_;

        Scalar tau_min_ = 0.99;
        Index max_iter_ = 200;
        Scalar tol_ = 1e-8;
        Scalar constr_viol_tol_ = 1e-8;
        Scalar mu_min_ = 1e-12;
        Scalar sigma_max_ = 1.0;
        bool verbose_ = true;

        Index iteration_ = 0;
        Scalar mu_ = 1.0;
        Scalar last_mu_aff_ = 0.;
    };

    /**
     * @brief Builder for MehrotraQpAlgorithm.
     *
     * Same role as @c IpAlgBuilder: takes an @c Nlp<ProblemType> (typically
     * @c NlpOcp wrapping an @c OcpAbstract) and wires up the IpData / NLP
     * wrapper / primal-dual solver / initializer stack used internally.
     */
    template <typename ProblemType> class MehrotraQpBuilder
    {
    public:
        explicit MehrotraQpBuilder(const std::shared_ptr<Nlp<ProblemType>> &nlp);

        MehrotraQpBuilder &with_options_registry(OptionRegistry *options_registry)
        {
            options_registry_ = options_registry;
            return *this;
        }

        std::shared_ptr<MehrotraQpAlgorithm<ProblemType>> build();

        std::shared_ptr<IpData<ProblemType>> get_ipdata() { return ipdata_; }

    private:
        std::shared_ptr<IpNlpOrig<ProblemType>> nlp_orig_;
        std::shared_ptr<IpData<ProblemType>> ipdata_;
        std::shared_ptr<ProblemInfo<ProblemType>> problem_info_;
        std::shared_ptr<AugSystemSolver<ProblemType>> aug_system_solver_;
        std::shared_ptr<PdSolverOrig<ProblemType>> pd_solver_;
        std::shared_ptr<IpInitializer<ProblemType>> initializer_;
        std::shared_ptr<IpEqMultInitializer<ProblemType>> eq_mult_initializer_;
        std::shared_ptr<MehrotraQpAlgorithm<ProblemType>> algorithm_;
        OptionRegistry *options_registry_ = nullptr;
    };

} // namespace fatrop

#endif // __fatrop_qp_fatrop_qp_hpp__
