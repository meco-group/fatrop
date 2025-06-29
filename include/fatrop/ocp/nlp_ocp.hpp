//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_nlp_ocp_hpp__
#define __fatrop_ocp_nlp_ocp_hpp__

#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/nlp/dims.hpp"
#include "fatrop/nlp/fwd.hpp"
#include "fatrop/nlp/nlp.hpp"
#include "fatrop/ocp/fwd.hpp"
#include "fatrop/ocp/ocp_abstract.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include "fatrop/ocp/type.hpp"
#include <memory>

namespace fatrop
{
    // OcpEvalType can be Eval abstract (dynamic polymorphism) or
    // a static implementation, a template specialization on the NlpOcpTpl level can be made as
    // well, using a newly created tag class. 
    template <typename OcpAbstractTag> class NlpOcpTpl : public Nlp<OcpType>
    {
        typedef std::shared_ptr<OcpAbstractTpl<OcpAbstractTag>> OcpAbstractSp;
        typedef ProblemInfo<OcpType> OcpInfo;

    public:
        NlpOcpTpl(const OcpAbstractSp &ocp);
        const NlpDims &nlp_dims() const override { return nlp_dims_; }
        const ProblemDims<OcpType> &problem_dims() const override { return ocp_dims_; }
        Index eval_lag_hess(const OcpInfo &info, const Scalar objective_scale,
                            const VecRealView &primal_x, const VecRealView &primal_s,
                            const VecRealView &lam, Hessian<OcpType> &hess) override;
        Index eval_constr_jac(const OcpInfo &info, const VecRealView &primal_x,
                              const VecRealView &primal_s, Jacobian<OcpType> &jac) override;
        Index eval_constraint_violation(const OcpInfo &info, const VecRealView &primal_x,
                                        const VecRealView &primal_s, VecRealView &res) override;
        Index eval_objective_gradient(const OcpInfo &info, const Scalar objective_scale,
                                      const VecRealView &primal_x, const VecRealView &primal_s,
                                      VecRealView &grad_x, VecRealView &grad_s) override;
        Index eval_objective(const OcpInfo &info, const Scalar objective_scale,
                             const VecRealView &primal_x, const VecRealView &primal_s,
                             Scalar &res) override;
        Index get_bounds(const OcpInfo &info, VecRealView &lower_bounds,
                         VecRealView &upper_bounds) override;
        Index get_initial_primal(const ProblemInfo<OcpType> &info, VecRealView &primal_x) override;
        void get_primal_damping(const ProblemInfo<OcpType> &info, VecRealView &damping) override;
        void apply_jacobian_s_transpose(const ProblemInfo<OcpType> &info,
                                                const VecRealView &multipliers, const Scalar alpha,
                                                const VecRealView &y, VecRealView &out) override;

    private:
        const OcpAbstractSp ocp_;
        ProblemDims<OcpType> ocp_dims_;
        NlpDims nlp_dims_;
    };

    typedef NlpOcpTpl<OcpAbstractDynamic> NlpOcp;

} // namespace fatrop

#endif //__fatrop_nlp_nlp__
