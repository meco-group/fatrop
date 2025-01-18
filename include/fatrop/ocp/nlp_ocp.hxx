//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_nlp_ocp_hxx__
#define __fatrop_ocp_nlp_ocp_hxx__

#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/nlp_ocp.hpp"
#include "fatrop/ocp/problem_info.hpp"

namespace fatrop
{
    namespace internal
    {
        template <typename OcpAbstractTag> struct NlpOcpAuxiliary
        {
            static OcpDims get_ocp_dims(const OcpAbstractTpl<OcpAbstractTag> &ocp)
            {
                const Index K = ocp.get_horizon_length();
                std::vector<Index> nu(K), nx(K), ng(K), ng_ineq(K);
                for (Index k = 0; k < K; k++)
                {
                    nu[k] = ocp.get_nu(k);
                    nx[k] = ocp.get_nx(k);
                    ng[k] = ocp.get_ng(k);
                    ng_ineq[k] = ocp.get_ng_ineq(k);
                }
                return OcpDims(K, nu, nx, ng, ng_ineq);
            }
            static NlpDims get_nlp_dims(const OcpDims &ocp_dims)
            {
                Index number_of_variables = 0;
                Index number_of_eq_constraints = 0;
                Index number_of_ineq_constraints = 0;
                for (Index k = 0; k < ocp_dims.K; k++)
                {
                    number_of_variables +=
                        ocp_dims.number_of_controls[k] + ocp_dims.number_of_states[k];
                    number_of_eq_constraints += ocp_dims.number_of_eq_constraints[k] +
                                                ocp_dims.number_of_ineq_constraints[k];
                    if (k != ocp_dims.K - 1)
                    {
                        number_of_eq_constraints += ocp_dims.number_of_states[k + 1];
                    }
                    number_of_ineq_constraints += ocp_dims.number_of_ineq_constraints[k];
                }
                return NlpDims(number_of_variables, number_of_eq_constraints,
                               number_of_ineq_constraints);
            }
        };
    }

    template <typename OcpAbstractTag>
    NlpOcpTpl<OcpAbstractTag>::NlpOcpTpl(const OcpAbstractSp &ocp)
        : ocp_(ocp),
          ocp_dims_(fatrop::internal::NlpOcpAuxiliary<OcpAbstractTag>::get_ocp_dims(*ocp)),
          nlp_dims_(fatrop::internal::NlpOcpAuxiliary<OcpAbstractTag>::get_nlp_dims(ocp_dims_))
    {
    }

    template <typename OcpAbstractTag>
    Index NlpOcpTpl<OcpAbstractTag>::eval_lag_hess(const OcpInfo &info,
                                                   const Scalar objective_scale,
                                                   const VecRealView &primal_x,
                                                   const VecRealView &primal_s,
                                                   const VecRealView &lam, Hessian<OcpType> &hess)
    {
        const Scalar * primal_x_p = primal_x.data();
        const Scalar * primal_s_p = primal_s.data();
        const Scalar * lam_p = lam.data();
        for (Index k = 0; k < info.dims.K; k++)
        {
            const Scalar *inputs_k = primal_x_p + info.offsets_primal_u[k];
            const Scalar *states_k = primal_x_p + info.offsets_primal_x[k];
            const Scalar *lam_dyn_k = lam_p + info.offsets_g_eq_dyn[k];
            const Scalar *lam_eq_k = lam_p + info.offsets_g_eq_path[k];
            const Scalar *lam_eq_ineq_k = lam_p + info.offsets_g_eq_slack[k];
            ocp_->eval_RSQrqt(&objective_scale, inputs_k, states_k, lam_dyn_k, lam_eq_k,
                              lam_eq_ineq_k, &hess.RSQrqt[k].mat(), k);
        }
        return 0;
    }

    template <typename OcpAbstractTag>
    Index
    NlpOcpTpl<OcpAbstractTag>::eval_constr_jac(const OcpInfo &info, const VecRealView &primal_x,
                                               const VecRealView &primal_s, Jacobian<OcpType> &jac)
    {
        const Scalar * primal_x_p = primal_x.data();
        for (Index k = 0; k < info.dims.K; k++)
        {
            const Scalar *inputs_k = primal_x_p + info.offsets_primal_u[k];
            const Scalar *states_k = primal_x_p + info.offsets_primal_x[k];
            ocp_->eval_Ggt(inputs_k, states_k, &jac.Gg_eqt[k].mat(), k);
            ocp_->eval_Ggt_ineq(inputs_k, states_k, &jac.Gg_ineqt[k].mat(), k);
            if (k != info.dims.K - 1)
            {
                const Scalar *states_kp1 = primal_x_p + info.offsets_primal_x[k + 1];
                ocp_->eval_BAbt(states_kp1, inputs_k, states_k, &jac.BAbt[k].mat(), k);
            }
        }
        return 0;
    }

    template <typename OcpAbstractTag>
    Index NlpOcpTpl<OcpAbstractTag>::eval_constraint_violation(const OcpInfo &info,
                                                               const VecRealView &primal_x,
                                                               const VecRealView &primal_s,
                                                               VecRealView &res)
    {
        const Scalar * primal_x_p = primal_x.data();
        Scalar * res_p = res.data();
        for (Index k = 0; k < info.dims.K; k++)
        {
            const Scalar *inputs_k = primal_x_p + info.offsets_primal_u[k];
            const Scalar *states_k = primal_x_p + info.offsets_primal_x[k];
            ocp_->eval_g(inputs_k, states_k, res_p + info.offsets_g_eq_path[k], k);
            ocp_->eval_gineq(inputs_k, states_k, res_p + info.offsets_g_eq_slack[k], k);
            if (k != info.dims.K - 1)
            {
                const Scalar *states_kp1 = primal_x_p + info.offsets_primal_x[k + 1];
                ocp_->eval_b(states_kp1, inputs_k, states_k, res_p + info.offsets_g_eq_dyn[k], k);
            }
        }
        // add -s to the slack constraints
        res.block(info.number_of_g_eq_slack, info.offset_g_eq_slack) =
            res.block(info.number_of_g_eq_slack, info.offset_g_eq_slack) -
            primal_s.block(info.number_of_g_eq_slack, 0);
        return 0;
    }

    template <typename OcpAbstractTag>
    Index NlpOcpTpl<OcpAbstractTag>::eval_objective_gradient(const OcpInfo &info,
                                                             const Scalar objective_scale,
                                                             const VecRealView &primal_x,
                                                             VecRealView &grad_x,
                                                             VecRealView &grad_s)
    {
        grad_s.block(info.number_of_g_eq_slack, 0) = 0;
        const Scalar * primal_x_p = primal_x.data();
        Scalar * grad_x_p = grad_x.data();
        for (Index k = 0; k < info.dims.K; k++)
        {
            const Scalar *inputs_k = primal_x_p + info.offsets_primal_u[k];
            const Scalar *states_k = primal_x_p + info.offsets_primal_x[k];
            ocp_->eval_rq(&objective_scale, inputs_k, states_k, grad_x_p + info.offsets_primal_u[k],
                          k);
        }
        return 0;
    }

    template <typename OcpAbstractTag>
    Index NlpOcpTpl<OcpAbstractTag>::eval_objective(const OcpInfo &info,
                                                    const Scalar objective_scale,
                                                    const VecRealView &primal_x,
                                                    const VecRealView &primal_s, Scalar &res)
    {
        res = 0;
        const Scalar * primal_x_p = primal_x.data();
        for (Index k = 0; k < info.dims.K; k++)
        {
            Scalar ret = 0;
            const Scalar *inputs_k = primal_x_p + info.offsets_primal_u[k];
            const Scalar *states_k = primal_x_p + info.offsets_primal_x[k];
            ocp_->eval_L(&objective_scale, inputs_k, states_k, &ret, k);
            res += ret;
        }
        return 0;
    }
}

#endif //__fatrop_ocp_nlp_ocp_hxx__
