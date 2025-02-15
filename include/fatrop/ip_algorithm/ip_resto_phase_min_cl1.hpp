//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_resto_phase_min_cl1_hpp__
#define __fatrop_ip_resto_phase_min_cl1_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/ip_algorithm/ip_resto_phase.hpp"
#include <memory>

namespace fatrop
{
    template <typename ProblemType> class IpRestoPhaseMinCl1 : public IpRestoPhaseBase
    {
        typedef IpAlgorithm<ProblemType> IpAlgorithmType;
        typedef std::shared_ptr<IpAlgorithmType> IpAlgorithmSp;
        typedef std::shared_ptr<IpEqMultInitializer<ProblemType>> IpEqMultInitializerSp;
        typedef IpNlpResto<ProblemType> IpNlpRestoType;
        typedef std::shared_ptr<IpNlpRestoType> IpNlpRestoSp;
        typedef IpIterate<ProblemType> IpIterateType;
        typedef IpData<ProblemType> IpDataType;
        typedef std::shared_ptr<IpDataType> IpDataSp;
        typedef ProblemInfo<ProblemType> ProblemInfoType;

    public:
        IpRestoPhaseMinCl1(const IpAlgorithmSp &ip_algorithm,
                           const IpEqMultInitializerSp &eq_mult_initializer,
                           const IpDataSp &ip_data_orig, const IpDataSp &ip_data_resto,
                           const IpNlpRestoSp &ip_nlp_resto);
        bool perform_restoration() override;
        void register_options(OptionRegistry &registry) override;
        // the linesearch is responsible for resetting the restoration phase
        void reset() override;

        virtual ~IpRestoPhaseMinCl1() = default;

        // Setters for options
        void set_bound_mult_reset_treshold(const Scalar &value)
        {
            bound_mult_reset_treshold_ = value;
        }
        void set_constr_mult_reset_treshold(const Scalar &value)
        {
            constr_mult_reset_treshold_ = value;
        }
        void set_resto_failure_feasibility_treshold(const Scalar &value)
        {
            resto_failure_feasibility_treshold_ = value;
        }
        void set_constr_viol_tol(const Scalar &value) { constr_viol_tol_ = value; }

    private:
        IpAlgorithmSp resto_ip_algorithm_;
        IpEqMultInitializerSp eq_mult_initializer_;
        IpDataSp ip_data_orig_;
        IpDataSp ip_data_resto_;
        IpNlpRestoSp ip_nlp_resto_;
        // options
        Scalar bound_mult_reset_treshold_ = 1e3;
        Scalar constr_mult_reset_treshold_ = 0.;
        Scalar resto_failure_feasibility_treshold_ = 0.;
        Scalar constr_viol_tol_ = 1e-4;
        // internal statistics
        Index count_resto_ = 0;
    };
} // namespace fatrop

#endif // __fatrop_ip_resto_phase_hpp__
