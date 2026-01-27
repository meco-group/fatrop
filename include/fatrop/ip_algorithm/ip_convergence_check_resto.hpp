//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_convergence_check_resto_hpp__
#define __fatrop_ip_convergence_check_resto_hpp__

#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/nlp/fwd.hpp"
#include "fatrop/ip_algorithm/ip_convergence_check.hpp"
namespace fatrop
{
    template <typename ProblemType>
    class IpConvergenceCheckResto : public IpConvergenceCheck<ProblemType>
    {
        typedef IpConvergenceCheck<ProblemType> Base;
        typedef IpData<ProblemType> IpDataType;
        typedef std::shared_ptr<IpDataType> IpDataSp;
        typedef ProblemInfo<ProblemType> ProblemInfoType;
        typedef IpIterate<ProblemType> IpIterateType;
        typedef std::shared_ptr<IpLineSearchBase> IpLineSearchSp;
        typedef std::weak_ptr<IpLineSearchBase> IpLineSearchWp;

    public:
        IpConvergenceCheckResto(const IpDataSp &data_orig, const IpDataSp &data_resto);
        IpConvergenceStatus check_converged() override;
        void register_options(OptionRegistry &registry) override;
        void reset() override;
        inline void set_required_infeasibility_reduction(const Scalar &val);
        inline void set_constr_viol_tol(const Scalar &val);
        inline void set_line_search_orig(const IpLineSearchSp &line_search);

    private:
        Scalar kappa_resto_ = 0.9;
        Scalar constr_viol_tol_ = 1e-4;
        bool first_resto_iter_;
        Index successive_resto_iter_;
        IpDataSp data_orig_;
        IpDataSp data_resto_;
        IpLineSearchWp line_search_orig_;
    };

    template <typename ProblemType>
    void
    IpConvergenceCheckResto<ProblemType>::set_required_infeasibility_reduction(const Scalar &val)
    {
        kappa_resto_ = val;
    }
    template <typename ProblemType>
    void IpConvergenceCheckResto<ProblemType>::set_constr_viol_tol(const Scalar &val)
    {
        constr_viol_tol_ = val;
    }
    template <typename ProblemType>
    void IpConvergenceCheckResto<ProblemType>::set_line_search_orig(const IpLineSearchSp &line_search)
    {
        line_search_orig_ = line_search;
    }

} // namespace fatrop

#endif // __fatrop_ip_convergence_check_resto_hpp__