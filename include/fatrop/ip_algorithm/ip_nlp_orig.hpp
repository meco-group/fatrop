//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_nlp_orig_hpp__
#define __fatrop_ip_algorithm_ip_nlp_orig_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include "fatrop/nlp/fwd.hpp"
#include "fatrop/nlp/nlp.hpp"
#include <memory>
namespace fatrop
{
    template <typename ProblemType> class IpNlpOrig : public Nlp<ProblemType>
    {
        typedef std::shared_ptr<Nlp<ProblemType>> NlpSp;

    public:
        IpNlpOrig(const NlpSp &nlp);
        virtual const NlpDims &nlp_dims() const override;
        virtual const ProblemDims<ProblemType> &problem_dims() const override;
        virtual Index eval_lag_hess(const ProblemInfo<ProblemType> &info,
                                    const Scalar objective_scale, const VecRealView &primal_x,
                                    const VecRealView &primal_s, const VecRealView &lam,
                                    Hessian<ProblemType> &hess) override;
        virtual Index eval_constr_jac(const ProblemInfo<ProblemType> &info,
                                      const VecRealView &primal_x, const VecRealView &primal_s,
                                      Jacobian<ProblemType> &jac) override;
        virtual Index eval_constraint_violation(const ProblemInfo<ProblemType> &info,
                                                const VecRealView &primal_x,
                                                const VecRealView &primal_s,
                                                VecRealView &res) override;
        virtual Index eval_objective_gradient(const ProblemInfo<ProblemType> &info,
                                              const Scalar objective_scale,
                                              const VecRealView &primal_x, VecRealView &grad_x,
                                              VecRealView &grad_s) override;
        virtual Index eval_objective(const ProblemInfo<ProblemType> &info,
                                     const Scalar objective_scale, const VecRealView &primal_x,
                                     const VecRealView &primal_s, Scalar &res) override;
        virtual Index get_bounds(const ProblemInfo<ProblemType> &info, VecRealView &lower_bounds,
                                 VecRealView &upper_bounds) override;

    private:
        NlpSp nlp_;
        VecRealAllocated modified_bounds_lower_;
        VecRealAllocated modified_bounds_upper_;
        Scalar constr_viol_tol_ = 1e-4;
        Scalar bound_relax_factor_ = 1e-8;
    };

} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_nlp_orig_hpp__
