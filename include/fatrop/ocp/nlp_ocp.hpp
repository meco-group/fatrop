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
        const OcpDims &ocp_dims() const { return ocp_dims_; }
        Index eval_lag_hess(const OcpInfo &info, const Scalar objective_scale,
                            const VecRealView &primal_x, const VecRealView &primal_s,
                            const VecRealView &lam, Hessian<OcpType> &hess) override;
        Index eval_constr_jac(const OcpInfo &info, const VecRealView &primal_x,
                              const VecRealView &primal_s, Jacobian<OcpType> &jac) override;
        Index eval_constraint_violation(const OcpInfo &info, const VecRealView &primal_x,
                                        const VecRealView &primal_s, VecRealView &res) override;
        Index eval_objective_gradient(const OcpInfo &info, const Scalar objective_scale,
                                      const VecRealView &primal_x, VecRealView &grad_x,
                                      VecRealView &grad_s) override;
        Index eval_objective(const OcpInfo &info, const Scalar objective_scale,
                             const VecRealView &primal_x, const VecRealView &primal_s,
                             Scalar &res) override;

    private:
        const OcpAbstractSp ocp_;
        OcpDims ocp_dims_;
        NlpDims nlp_dims_;
    };

    typedef NlpOcpTpl<OcpAbstractDynamic> NlpOcp;

} // namespace fatrop

#endif //__fatrop_nlp_nlp__