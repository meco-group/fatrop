//
// Copyright (c) Lander Vanroye, KU Leuven
//
#include "fatrop/ocp/jacobian.hpp"
#include "fatrop/common/exception.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include "fatrop/ocp/dims.hpp"
#include "fatrop/ocp/problem_info.hpp"
using namespace fatrop;

Jacobian<OcpType>::Jacobian(const OcpDims &dims)
{
    // reserve memory for the Jacobian matrices
    BAbt.reserve(dims.K - 1);
    Gg_eqt.reserve(dims.K);
    Gg_ineqt.reserve(dims.K);
    // allocate memory for the Jacobian matrices
    for (int k = 0; k < dims.K - 1; ++k)
        BAbt.emplace_back(dims.number_of_states[k] + dims.number_of_controls[k] + 1,
                          dims.number_of_states[k + 1]);
    for (int k = 0; k < dims.K; ++k)
    {
        Gg_eqt.emplace_back(dims.number_of_states[k] + dims.number_of_controls[k] + 1,
                            dims.number_of_eq_constraints[k]);
    }
    for (int k = 0; k < dims.K; ++k)
        Gg_ineqt.emplace_back(dims.number_of_states[k] + dims.number_of_controls[k] + 1,
                              dims.number_of_ineq_constraints[k]);
};
void Jacobian<OcpType>::apply_on_right(const OcpInfo &info, const VecRealView &x, Scalar alpha,
                                       const VecRealView &y, VecRealView &out) const
{
    out = alpha*y;
    // dynamics constraints BA*ux - x_next
    for (Index k = 0; k < info.dims.K - 1; ++k)
    {
        Index nu = info.dims.number_of_controls[k];
        Index nx = info.dims.number_of_states[k];
        Index nx_next = info.dims.number_of_states[k + 1];
        Index offset_ux = info.offsets_primal_u[k];
        Index offset_x_next = info.offsets_primal_x[k + 1];
        Index offset_dyn_eq = info.offsets_g_eq_dyn[k];
        // apply out[offs:offs+nx] =  BAbt.T @ x[offs:offs+nu+nx] - x_next[offs:offs+nx]
        gemv_t(nu + nx, nx_next, 1.0, BAbt[k], 0, 0, x, offset_ux, 1.0, out, offset_dyn_eq, out,
               offset_dyn_eq);
        axpy(nx_next, -1.0, x, offset_x_next, out, offset_dyn_eq, out, offset_dyn_eq);
    }
    // equality path constraints
    for (Index k = 0; k < info.dims.K; ++k)
    {
        Index nu = info.dims.number_of_controls[k];
        Index nx = info.dims.number_of_states[k];
        Index ng = info.dims.number_of_eq_constraints[k];
        Index offset_ux = info.offsets_primal_u[k];
        Index offset_g_eq = info.offsets_g_eq_path[k];
        // apply out[offs:offs+ng] =  Gg_eqt.T @ x[offs:offs+nu+nx]
        gemv_t(nu + nx, ng, 1.0, Gg_eqt[k], 0, 0, x, offset_ux, 1.0, out, offset_g_eq, out,
               offset_g_eq);
    }
    // slack equality path constraints Gg_ineqt @ x
    for (Index k = 0; k < info.dims.K; ++k)
    {
        Index nu = info.dims.number_of_controls[k];
        Index nx = info.dims.number_of_states[k];
        Index ng_ineq = info.dims.number_of_ineq_constraints[k];
        Index offset_ux = info.offsets_primal_u[k];
        Index offset_g_ineq = info.offsets_g_eq_slack[k];
        // apply out[offs:offs+ng_ineq] =  Gg_ineqt.T @ x[offs:offs+nu+nx]
        gemv_t(nu + nx, ng_ineq, 1.0, Gg_ineqt[k], 0, 0, x, offset_ux, 1.0, out, offset_g_ineq, out,
               offset_g_ineq);
    }
};
void Jacobian<OcpType>::transpose_apply_on_right(const OcpInfo &info, const VecRealView &mult_eq,
                                                 Scalar alpha, const VecRealView &y,
                                                 VecRealView &out) const
{
    // set the output to zero, we will add the contributions
    // dynamics constraints'contributions [BA.T ; 0; -I] @ mult_eq
    out = alpha*y;
    for (Index k = 0; k < info.dims.K - 1; ++k)
    {
        Index nu = info.dims.number_of_controls[k];
        Index nx = info.dims.number_of_states[k];
        Index offs_ux = info.offsets_primal_u[k];
        Index offs_x_next = info.offsets_primal_x[k + 1];
        Index nx_next = info.dims.number_of_states[k + 1];
        Index offset_g_dyn = info.offsets_g_eq_dyn[k];
        // apply out[offs_ux:offs_ux + nu + nx] +=  BAbt @ mult_eq[offs_g_dyn:offs_g_dyn + nx_next]
        gemv_n(nu + nx, nx_next, 1.0, BAbt[k], 0, 0, mult_eq, offset_g_dyn, 1.0, out, offs_ux, out,
               offs_ux);
        // apply out[offs_ux:offs_ux + nu + nx] -= mult_eq[offs_f_dyn:offs_g_dyn + nx_next]
        axpy(nx_next, -1.0, mult_eq, offset_g_dyn, out, offs_x_next, out, offs_x_next);
    };
    // equality path constraints' contributions Gg_eqt @ mult_eq
    for (Index k = 0; k < info.dims.K; ++k)
    {
        Index nu = info.dims.number_of_controls[k];
        Index nx = info.dims.number_of_states[k];
        Index ng = info.dims.number_of_eq_constraints[k];
        Index offset_ux = info.offsets_primal_u[k];
        Index offset_g_eq = info.offsets_g_eq_path[k];
        // apply out[offs:offs+nu+nx] +=  Gg_eqt @ mult_eq[offs:offs+ng]
        gemv_n(nu + nx, ng, 1.0, Gg_eqt[k], 0, 0, mult_eq, offset_g_eq, 1.0, out, offset_ux, out,
               offset_ux);
    }
    // inequality path constraints' contributions Gg_ineqt @ mult_eq
    for (Index k = 0; k < info.dims.K; ++k)
    {
        Index nu = info.dims.number_of_controls[k];
        Index nx = info.dims.number_of_states[k];
        Index ng_ineq = info.dims.number_of_ineq_constraints[k];
        Index offset_ux = info.offsets_primal_u[k];
        Index offset_g_ineq = info.offsets_g_eq_slack[k];
        // apply out[offs:offs+nu+nx] +=  Gg_ineqt @ mult_eq[offs:offs+ng_ineq]
        gemv_n(nu + nx, ng_ineq, 1.0, Gg_ineqt[k], 0, 0, mult_eq, offset_g_ineq, 1.0, out,
               offset_ux, out, offset_ux);
    }
}
void Jacobian<OcpType>::get_rhs(const OcpInfo &info, VecRealView &rhs) const
{
    // dynamics constraints' right-hand side
    for (Index k = 0; k < info.dims.K - 1; ++k)
    {
        Index nu = info.dims.number_of_controls[k];
        Index nx = info.dims.number_of_states[k];
        Index nx_next = info.dims.number_of_states[k + 1];
        Index offset_eq_dyn = info.offsets_g_eq_dyn[k];
        // the rhs is the last row of the BAbt[k] matrix
        rowex(nx_next, 1.0, BAbt[k], nu + nx, 0, rhs, offset_eq_dyn);
    }
    // equality path constraints' right-hand side
    for (Index k = 0; k < info.dims.K; ++k)
    {
        Index nu = info.dims.number_of_controls[k];
        Index nx = info.dims.number_of_states[k];
        Index ng = info.dims.number_of_eq_constraints[k];
        Index offset_eq_path = info.offsets_g_eq_path[k];
        // the rhs is the last row of the Gg_eqt[k] matrix
        rowex(ng, 1.0, Gg_eqt[k], nu + nx, 0, rhs, offset_eq_path);
    }
    // inequality path constraints' right-hand side
    for (Index k = 0; k < info.dims.K; ++k)
    {
        Index nu = info.dims.number_of_controls[k];
        Index nx = info.dims.number_of_states[k];
        Index ng_ineq = info.dims.number_of_ineq_constraints[k];
        Index offset_eq_ineq = info.offsets_g_eq_slack[k];
        // the rhs is the last row of the Gg_ineqt[k] matrix
        rowex(ng_ineq, 1.0, Gg_ineqt[k], nu + nx, 0, rhs, offset_eq_ineq);
    }
};
void Jacobian<OcpType>::set_rhs(const OcpInfo &info, const VecRealView &rhs)
{
    // dynamics constraints' right-hand side
    for (Index k = 0; k < info.dims.K - 1; ++k)
    {
        Index nu = info.dims.number_of_controls[k];
        Index nx = info.dims.number_of_states[k];
        Index nx_next = info.dims.number_of_states[k + 1];
        Index offset_eq_dyn = info.offsets_g_eq_dyn[k];
        // the rhs is the last row of the BAbt[k] matrix
        rowin(nx_next, 1.0, rhs, offset_eq_dyn, BAbt[k], nu + nx, 0);
    }
    // equality path constraints' right-hand side
    for (Index k = 0; k < info.dims.K; ++k)
    {
        Index nu = info.dims.number_of_controls[k];
        Index nx = info.dims.number_of_states[k];
        Index ng = info.dims.number_of_eq_constraints[k];
        Index offset_eq_path = info.offsets_g_eq_path[k];
        // the rhs is the last row of the Gg_eqt[k] matrix
        rowin(ng, 1.0, rhs, offset_eq_path, Gg_eqt[k], nu + nx, 0);
    }
    // inequality path constraints' right-hand side
    for (Index k = 0; k < info.dims.K; ++k)
    {
        Index nu = info.dims.number_of_controls[k];
        Index nx = info.dims.number_of_states[k];
        Index ng_ineq = info.dims.number_of_ineq_constraints[k];
        Index offset_eq_ineq = info.offsets_g_eq_slack[k];
        // the rhs is the last row of the Gg_ineqt[k] matrix
        rowin(ng_ineq, 1.0, rhs, offset_eq_ineq, Gg_ineqt[k], nu + nx, 0);
    }
};