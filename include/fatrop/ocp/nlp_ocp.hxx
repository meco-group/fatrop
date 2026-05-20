//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_nlp_ocp_hxx__
#define __fatrop_ocp_nlp_ocp_hxx__

#include "fatrop/ocp/hessian.hpp"
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/ocp/nlp_ocp.hpp"
#include "fatrop/ocp/problem_info.hpp"
#include <type_traits>

namespace fatrop
{
    namespace internal
    {
        // SFINAE detection for the tangent / retraction extensions added on top of the
        // OcpAbstractTpl interface. Static (compile-time) OCP specializations are not
        // required to override these methods — the helpers below fall back to the
        // Euclidean defaults (matching the dynamic base class).
        template <typename, typename = void> struct has_get_nx_tangent : std::false_type {};
        template <typename T>
        struct has_get_nx_tangent<
            T, std::void_t<decltype(std::declval<const T &>().get_nx_tangent(Index{}))>>
            : std::true_type {};
        template <typename, typename = void> struct has_get_nu_tangent : std::false_type {};
        template <typename T>
        struct has_get_nu_tangent<
            T, std::void_t<decltype(std::declval<const T &>().get_nu_tangent(Index{}))>>
            : std::true_type {};
        template <typename, typename = void> struct has_apply_retraction_xk : std::false_type {};
        template <typename T>
        struct has_apply_retraction_xk<
            T, std::void_t<decltype(std::declval<T &>().apply_retraction_xk(
                   Index{}, (const Scalar *)nullptr, (const Scalar *)nullptr, Scalar{},
                   (Scalar *)nullptr))>> : std::true_type {};
        template <typename, typename = void> struct has_apply_retraction_uk : std::false_type {};
        template <typename T>
        struct has_apply_retraction_uk<
            T, std::void_t<decltype(std::declval<T &>().apply_retraction_uk(
                   Index{}, (const Scalar *)nullptr, (const Scalar *)nullptr, Scalar{},
                   (Scalar *)nullptr))>> : std::true_type {};

        template <typename Ocp> Index call_get_nx_tangent(const Ocp &ocp, const Index k)
        {
            if constexpr (has_get_nx_tangent<Ocp>::value)
                return ocp.get_nx_tangent(k);
            else
                return ocp.get_nx(k);
        }
        template <typename Ocp> Index call_get_nu_tangent(const Ocp &ocp, const Index k)
        {
            if constexpr (has_get_nu_tangent<Ocp>::value)
                return ocp.get_nu_tangent(k);
            else
                return ocp.get_nu(k);
        }
        template <typename Ocp>
        void call_apply_retraction_xk(Ocp &ocp, const Index k, const Scalar *xk,
                                      const Scalar *delta_xk, const Scalar alpha,
                                      Scalar *xk_next)
        {
            if constexpr (has_apply_retraction_xk<Ocp>::value)
            {
                ocp.apply_retraction_xk(k, xk, delta_xk, alpha, xk_next);
            }
            else
            {
                const Index n = ocp.get_nx(k);
                for (Index i = 0; i < n; ++i)
                    xk_next[i] = xk[i] + alpha * delta_xk[i];
            }
        }
        template <typename Ocp>
        void call_apply_retraction_uk(Ocp &ocp, const Index k, const Scalar *uk,
                                      const Scalar *delta_uk, const Scalar alpha,
                                      Scalar *uk_next)
        {
            if constexpr (has_apply_retraction_uk<Ocp>::value)
            {
                ocp.apply_retraction_uk(k, uk, delta_uk, alpha, uk_next);
            }
            else
            {
                const Index n = ocp.get_nu(k);
                for (Index i = 0; i < n; ++i)
                    uk_next[i] = uk[i] + alpha * delta_uk[i];
            }
        }

        template <typename, typename = void>
        struct has_apply_dual_eq_dyn_transformation_k : std::false_type {};
        template <typename T>
        struct has_apply_dual_eq_dyn_transformation_k<
            T, std::void_t<decltype(std::declval<T &>().apply_dual_eq_dyn_transformation_k(
                   Index{}, (const Scalar *)nullptr, (const Scalar *)nullptr,
                   (const Scalar *)nullptr, (const Scalar *)nullptr,
                   (Scalar *)nullptr))>> : std::true_type {};
        template <typename, typename = void>
        struct has_apply_dual_eq_path_transformation_k : std::false_type {};
        template <typename T>
        struct has_apply_dual_eq_path_transformation_k<
            T, std::void_t<decltype(std::declval<T &>().apply_dual_eq_path_transformation_k(
                   Index{}, (const Scalar *)nullptr, (const Scalar *)nullptr,
                   (const Scalar *)nullptr, (Scalar *)nullptr))>> : std::true_type {};
        template <typename, typename = void>
        struct has_apply_dual_eq_slack_transformation_k : std::false_type {};
        template <typename T>
        struct has_apply_dual_eq_slack_transformation_k<
            T, std::void_t<decltype(std::declval<T &>().apply_dual_eq_slack_transformation_k(
                   Index{}, (const Scalar *)nullptr, (const Scalar *)nullptr,
                   (const Scalar *)nullptr, (Scalar *)nullptr))>> : std::true_type {};

        template <typename Ocp>
        void call_apply_dual_eq_dyn_transformation_k(Ocp &ocp, const Index k, const Scalar *xk,
                                                    const Scalar *xkp1, const Scalar *uk,
                                                    const Scalar *dual_in, Scalar *dual_out)
        {
            if constexpr (has_apply_dual_eq_dyn_transformation_k<Ocp>::value)
            {
                ocp.apply_dual_eq_dyn_transformation_k(k, xk, xkp1, uk, dual_in, dual_out);
            }
            else
            {
                const Index n = call_get_nx_tangent(ocp, k + 1);
                for (Index i = 0; i < n; ++i)
                    dual_out[i] = dual_in[i];
            }
        }

        template <typename Ocp>
        void call_apply_dual_eq_path_transformation_k(Ocp &ocp, const Index k, const Scalar *uk,
                                                     const Scalar *xk, const Scalar *dual_in,
                                                     Scalar *dual_out)
        {
            if constexpr (has_apply_dual_eq_path_transformation_k<Ocp>::value)
            {
                ocp.apply_dual_eq_path_transformation_k(k, uk, xk, dual_in, dual_out);
            }
            else
            {
                const Index n = ocp.get_ng(k);
                for (Index i = 0; i < n; ++i)
                    dual_out[i] = dual_in[i];
            }
        }

        template <typename Ocp>
        void call_apply_dual_eq_slack_transformation_k(Ocp &ocp, const Index k, const Scalar *uk,
                                                      const Scalar *xk, const Scalar *dual_in,
                                                      Scalar *dual_out)
        {
            if constexpr (has_apply_dual_eq_slack_transformation_k<Ocp>::value)
            {
                ocp.apply_dual_eq_slack_transformation_k(k, uk, xk, dual_in, dual_out);
            }
            else
            {
                const Index n = ocp.get_ng_ineq(k);
                for (Index i = 0; i < n; ++i)
                    dual_out[i] = dual_in[i];
            }
        }

        template <typename OcpAbstractTag> struct NlpOcpAuxiliary
        {
            static ProblemDims<OcpType> get_ocp_dims(const OcpAbstractTpl<OcpAbstractTag> &ocp)
            {
                const Index K = ocp.get_horizon_length();
                std::vector<Index> nu(K), nx(K), ng(K), ng_ineq(K);
                std::vector<Index> nu_tan(K), nx_tan(K);
                for (Index k = 0; k < K; k++)
                {
                    nu[k] = ocp.get_nu(k);
                    nx[k] = ocp.get_nx(k);
                    ng[k] = ocp.get_ng(k);
                    ng_ineq[k] = ocp.get_ng_ineq(k);
                    nu_tan[k] = call_get_nu_tangent(ocp, k);
                    nx_tan[k] = call_get_nx_tangent(ocp, k);
                }
                return ProblemDims<OcpType>(K, nu, nx, ng, ng_ineq, nu_tan, nx_tan);
            }
            static NlpDims get_nlp_dims(const ProblemDims<OcpType> &ocp_dims)
            {
                Index number_of_variables = 0;
                Index number_of_tangent_variables = 0;
                Index number_of_eq_constraints = 0;
                Index number_of_ineq_constraints = 0;
                for (Index k = 0; k < ocp_dims.K; k++)
                {
                    number_of_variables +=
                        ocp_dims.number_of_controls[k] + ocp_dims.number_of_states[k];
                    number_of_tangent_variables += ocp_dims.number_of_tangent_controls[k] +
                                                   ocp_dims.number_of_tangent_states[k];
                    number_of_eq_constraints += ocp_dims.number_of_eq_constraints[k] +
                                                ocp_dims.number_of_ineq_constraints[k];
                    if (k != ocp_dims.K - 1)
                    {
                        // Dynamics constraint x_{k+1} = f(x_k, u_k) is linearized in tangent space,
                        // so the equality block contributes the tangent state size of stage k+1.
                        number_of_eq_constraints += ocp_dims.number_of_tangent_states[k + 1];
                    }
                    number_of_ineq_constraints += ocp_dims.number_of_ineq_constraints[k];
                }
                return NlpDims(number_of_variables, number_of_tangent_variables,
                               number_of_eq_constraints, number_of_ineq_constraints);
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
        const Scalar *primal_x_p = primal_x.data();
        const Scalar *primal_s_p = primal_s.data();
        const Scalar *lam_p = lam.data();
        for (Index k = 0; k < info.dims.K; k++)
        {
            const Scalar *inputs_k = primal_x_p + info.offsets_primal_u[k];
            const Scalar *states_k = primal_x_p + info.offsets_primal_x[k];
            const Scalar *lam_dyn_k =
                (k != info.dims.K - 1) ? lam_p + info.offsets_g_eq_dyn[k] : nullptr;
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
        const Scalar *primal_x_p = primal_x.data();
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
        const Scalar *primal_x_p = primal_x.data();
        Scalar *res_p = res.data();
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
    Index NlpOcpTpl<OcpAbstractTag>::eval_objective_gradient(
        const OcpInfo &info, const Scalar objective_scale, const VecRealView &primal_x,
        const VecRealView &primal_s, VecRealView &grad_x, VecRealView &grad_s)
    {
        grad_s.block(info.number_of_g_eq_slack, 0) = 0;
        const Scalar *primal_x_p = primal_x.data();
        Scalar *grad_x_p = grad_x.data();
        for (Index k = 0; k < info.dims.K; k++)
        {
            // Read state/control values from the *primal* (manifold) variable vector.
            const Scalar *inputs_k = primal_x_p + info.offsets_primal_u[k];
            const Scalar *states_k = primal_x_p + info.offsets_primal_x[k];
            // Write the gradient into the *tangent* (search-direction) vector. For
            // Euclidean problems the tangent offsets coincide with the primal offsets.
            ocp_->eval_rq(&objective_scale, inputs_k, states_k,
                          grad_x_p + info.offsets_tangent_u[k], k);
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
        const Scalar *primal_x_p = primal_x.data();
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
    template <typename OcpAbstractTag>
    Index NlpOcpTpl<OcpAbstractTag>::get_bounds(const OcpInfo &info, VecRealView &lower_bounds,
                                                VecRealView &upper_bounds)
    {
        if (info.number_of_slack_variables == 0)
            return 0;
        Scalar *lower_bounds_p = lower_bounds.data();
        Scalar *upper_bounds_p = upper_bounds.data();
        for (Index k = 0; k < info.dims.K; k++)
        {
            Scalar *lower_bounds_k = lower_bounds_p + info.offsets_slack[k];
            Scalar *upper_bounds_k = upper_bounds_p + info.offsets_slack[k];
            ocp_->get_bounds(lower_bounds_k, upper_bounds_k, k);
        }
        return 0;
    }

    template <typename OcpAbstractTag>
    Index NlpOcpTpl<OcpAbstractTag>::get_initial_primal(const ProblemInfo<OcpType> &info,
                                                        VecRealView &primal_x)
    {
        Scalar *primal_x_ptr = primal_x.data();
        for (Index k = 0; k < info.dims.K; k++)
        {
            ocp_->get_initial_uk(primal_x_ptr + info.offsets_primal_u[k], k);
            ocp_->get_initial_xk(primal_x_ptr + info.offsets_primal_x[k], k);
        }
        return 0;
    }
    template <typename OcpAbstractTag>
    void NlpOcpTpl<OcpAbstractTag>::get_primal_damping(const ProblemInfo<OcpType> &info,
                                                       VecRealView &damping)
    {
        damping = 0;
    }
    template <typename OcpAbstractTag>
    void NlpOcpTpl<OcpAbstractTag>::apply_jacobian_s_transpose(const ProblemInfo<OcpType> &info,
                                                               const VecRealView &multipliers,
                                                               const Scalar alpha,
                                                               const VecRealView &y,
                                                               VecRealView &out)
    {
        out = alpha * y;
        out.block(info.number_of_slack_variables, 0) =
            out.block(info.number_of_slack_variables, 0) -
            multipliers.block(info.number_of_slack_variables, info.offset_g_eq_slack);
    }

    template <typename OcpAbstractTag>
    void NlpOcpTpl<OcpAbstractTag>::apply_dual_eq_transformation(
        const ProblemInfo<OcpType> &info, const VecRealView &primal_x,
        const VecRealView &dual_eq_in, VecRealView &dual_eq_out)
    {
        const Scalar *primal_x_p = primal_x.data();
        const Scalar *dual_in_p = dual_eq_in.data();
        Scalar *dual_out_p = dual_eq_out.data();
        // Dynamics multipliers (one block per inter-stage transition).
        for (Index k = 0; k < info.dims.K - 1; k++)
        {
            const Scalar *xk = primal_x_p + info.offsets_primal_x[k];
            const Scalar *xkp1 = primal_x_p + info.offsets_primal_x[k + 1];
            const Scalar *uk = primal_x_p + info.offsets_primal_u[k];
            const Scalar *lin = dual_in_p + info.offsets_g_eq_dyn[k];
            Scalar *lout = dual_out_p + info.offsets_g_eq_dyn[k];
            internal::call_apply_dual_eq_dyn_transformation_k(*ocp_, k, xk, xkp1, uk, lin, lout);
        }
        // Equality-path multipliers per stage.
        for (Index k = 0; k < info.dims.K; k++)
        {
            const Scalar *xk = primal_x_p + info.offsets_primal_x[k];
            const Scalar *uk = primal_x_p + info.offsets_primal_u[k];
            const Scalar *lin = dual_in_p + info.offsets_g_eq_path[k];
            Scalar *lout = dual_out_p + info.offsets_g_eq_path[k];
            internal::call_apply_dual_eq_path_transformation_k(*ocp_, k, uk, xk, lin, lout);
        }
        // Slack (inequality-as-equality) multipliers per stage.
        for (Index k = 0; k < info.dims.K; k++)
        {
            const Scalar *xk = primal_x_p + info.offsets_primal_x[k];
            const Scalar *uk = primal_x_p + info.offsets_primal_u[k];
            const Scalar *lin = dual_in_p + info.offsets_g_eq_slack[k];
            Scalar *lout = dual_out_p + info.offsets_g_eq_slack[k];
            internal::call_apply_dual_eq_slack_transformation_k(*ocp_, k, uk, xk, lin, lout);
        }
    }

    template <typename OcpAbstractTag>
    void NlpOcpTpl<OcpAbstractTag>::apply_retraction(const ProblemInfo<OcpType> &info,
                                                     const VecRealView &primal_x,
                                                     const VecRealView &delta_primal_x,
                                                     const Scalar alpha,
                                                     VecRealView &primal_x_next)
    {
        const Scalar *primal_x_p = primal_x.data();
        const Scalar *delta_x_p = delta_primal_x.data();
        Scalar *primal_x_next_p = primal_x_next.data();
        for (Index k = 0; k < info.dims.K; k++)
        {
            const Scalar *uk_in = primal_x_p + info.offsets_primal_u[k];
            const Scalar *xk_in = primal_x_p + info.offsets_primal_x[k];
            const Scalar *duk_in = delta_x_p + info.offsets_tangent_u[k];
            const Scalar *dxk_in = delta_x_p + info.offsets_tangent_x[k];
            Scalar *uk_out = primal_x_next_p + info.offsets_primal_u[k];
            Scalar *xk_out = primal_x_next_p + info.offsets_primal_x[k];
            internal::call_apply_retraction_uk(*ocp_, k, uk_in, duk_in, alpha, uk_out);
            internal::call_apply_retraction_xk(*ocp_, k, xk_in, dxk_in, alpha, xk_out);
        }
    }
}

#endif //__fatrop_ocp_nlp_ocp_hxx__
