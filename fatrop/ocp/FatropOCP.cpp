/*
 * Fatrop - A fast trajectory optimization solver
 * Copyright (C) 2022, 2023 Lander Vanroye, KU Leuven. All rights reserved.
 *
 * This file is part of Fatrop.
 *
 * Fatrop is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fatrop is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Fatrop.  If not, see <http://www.gnu.org/licenses/>. */
#include "ocp/FatropOCP.hpp"
using namespace fatrop;
using namespace std;
FatropOCP::FatropOCP(
    const shared_ptr<OCP> &ocp,
    const shared_ptr<OCPLinearSolver> &ls,
    const shared_ptr<OCPScalingMethod> &scaler, const shared_ptr<FatropOptions> &options, const shared_ptr<FatropPrinter> &printer) : ocp_(ocp), dims_(ocp_->get_ocp_dims()),
                                                                                                                                      nlpdims_({sum(dims_.nx + dims_.nu), sum(dims_.ng + dims_.ng_ineq + dims_.nx) - dims_.nx.get(0), sum(dims_.ng_ineq)}), ls_(ls), scaler_(scaler), options_(options), printer_(printer), ocpkktmemory_(dims_), s_memvec(nlpdims_.nineqs, 4), ux_memvec(nlpdims_.nvars, 1),
                                                                                                                                      sigma(s_memvec[0]),
                                                                                                                                      gradb(s_memvec[1]),
                                                                                                                                      s_dummy(s_memvec[2]),
                                                                                                                                      s_zero(s_memvec[3]),
                                                                                                                                      ux_dummy(ux_memvec[0]),
                                                                                                                                      rhs_rq(nlpdims_.nvars, 1),
                                                                                                                                      rhs_b(dims_.n_b_tot, 1),
                                                                                                                                      rhs_g(dims_.n_g_tot, 1),
                                                                                                                                      rhs_g_ineq(dims_.n_g_ineq_tot, 1),
                                                                                                                                      rhs_gradb(dims_.n_g_ineq_tot, 1),
                                                                                                                                      rhs_rq2(nlpdims_.nvars, 1),
                                                                                                                                      rhs_b2(dims_.n_b_tot, 1),
                                                                                                                                      rhs_g2(dims_.n_g_tot, 1),
                                                                                                                                      rhs_g_ineq2(dims_.n_g_ineq_tot, 1),
                                                                                                                                      rhs_gradb2(dims_.n_g_ineq_tot, 1),
                                                                                                                                      gradb_total_cache(dims_.n_g_ineq_tot, 1),
                                                                                                                                      sigma_total_cache(dims_.n_g_ineq_tot, 1),
                                                                                                                                      ux_test(nlpdims_.nvars, 1),
                                                                                                                                      lam_test(nlpdims_.neqs, 1),
                                                                                                                                      delta_s_test(nlpdims_.nineqs, 1)
{
    options_->register_option(BooleanOption("iterative_refinement_SOC", "Use iterative refinement for SOC", &it_ref, true));
}
int FatropOCP::eval_lag_hess(
    double obj_scale,
    const FatropVecBF &primal_vars,
    const FatropVecBF &lam)
{
    int res = ocp_->eval_lag_hess(
        &ocpkktmemory_,
        obj_scale,
        primal_vars,
        lam);
    return res;
};
int FatropOCP::eval_constr_jac(
    const FatropVecBF &primal_vars,
    const FatropVecBF &slack_vars)
{
    int res = ocp_->eval_constr_jac(
        &ocpkktmemory_,
        primal_vars,
        slack_vars);
    return res;
};
int FatropOCP::solve_pd_sys(
    const double inertia_correction_w,
    const double inertia_correction_c,
    const FatropVecBF &ux,
    const FatropVecBF &lam,
    const FatropVecBF &delta_s,
    const FatropVecBF &sigma_total,
    const FatropVecBF &gradb_total)
{
    // ls_ = RefCountPtr<OCPLinearSolver>(new Sparse_OCP(ocp_->GetOCPDims(), ocpkktmemory_));
    // save gradb_total and sigma_total
    gradb_total_cache[0].copy(gradb_total);
    sigma_total_cache[0].copy(sigma_total);
    inertia_correction_w_cache = inertia_correction_w;
    inertia_correction_c_cache = inertia_correction_c;

    return ls_->solve_pd_sys(
        &ocpkktmemory_,
        inertia_correction_w,
        inertia_correction_c,
        ux,
        lam,
        delta_s,
        sigma_total,
        gradb_total);
};

int FatropOCP::solve_soc_rhs(
    const FatropVecBF &ux,
    const FatropVecBF &lam,
    const FatropVecBF &delta_s,
    const FatropVecBF &constraint_violation)
{
    int min_it_ref = 0;
    if (inertia_correction_c_cache != 0.0)
        return -1;
    // bool it_ref = true;
    /// todo avoid retrieving unnecessary rhs'es
    ls_->get_rhs(
        &ocpkktmemory_,
        gradb_total_cache[0],
        rhs_rq[0],
        rhs_b[0],
        rhs_g[0],
        rhs_g_ineq[0],
        rhs_gradb[0]);
    // prepare rhs_b, rhs_g and rhg_g_ineq
    rhs_g[0].copy(constraint_violation.block(0, dims_.n_g_tot));
    rhs_b[0].copy(constraint_violation.block(dims_.n_g_tot, dims_.n_b_tot));
    rhs_g_ineq[0].copy(constraint_violation.block(dims_.n_b_tot + dims_.n_g_tot, dims_.n_g_ineq_tot));

    ls_->solve_rhs(
        &ocpkktmemory_,
        ux,
        lam,
        delta_s,
        sigma_total_cache[0],
        rhs_rq[0],
        rhs_b[0],
        rhs_g[0],
        rhs_g_ineq[0],
        rhs_gradb[0]);
    if (it_ref)
    {
        double err_curr = 0.0;
        // copy(ux, ux_test[0]);
        // copy(lam, lam_test[0]);
        // copy(delta_s, delta_s_test[0]);
        // GetRHS(
        //     OCP,
        //     gradb_total,
        //     rhs_rq2[0],
        //     rhs_b2[0],
        //     rhs_g2[0],
        //     rhs_g_ineq2[0],
        //     rhs_gradb2[0]);
        double max_norm = max(Linf(rhs_gradb[0]), max(Linf(rhs_g_ineq[0]), max(Linf(rhs_g[0]), max(Linf(rhs_rq[0]), Linf(rhs_b[0])))));
        max_norm = (max_norm == 0.0) ? 1.0 : max_norm;
        double error_prev = -1.0;
        for (int i = 0; i < 5; i++)
        {
            ls_->compute_pd_sys_times_vec(
                &ocpkktmemory_,
                inertia_correction_w_cache,
                0.0,
                ux,
                lam,
                delta_s,
                sigma_total_cache[0],
                rhs_rq2[0],
                rhs_b2[0],
                rhs_g2[0],
                rhs_g_ineq2[0],
                rhs_gradb2[0]);
            axpby(1.0, rhs_rq2[0], 1.0, rhs_rq[0], rhs_rq2[0]);
            axpby(1.0, rhs_b2[0], 1.0, rhs_b[0], rhs_b2[0]);
            axpby(1.0, rhs_g2[0], 1.0, rhs_g[0], rhs_g2[0]);
            axpby(1.0, rhs_g_ineq2[0], 1.0, rhs_g_ineq[0], rhs_g_ineq2[0]);
            axpby(1.0, rhs_gradb2[0], 1.0, rhs_gradb[0], rhs_gradb2[0]);

            // cout << "residu rq:  " << Linf(rhs_rq[0]) / max_norm << "  ";
            // cout << "residu b:  " << Linf(rhs_b[0]) / max_norm << "  ";
            // cout << "residu g:  " << Linf(rhs_g[0]) / max_norm << "  ";
            // cout << "residu g_ineq:  " << Linf(rhs_g_ineq[0]) / max_norm << "  ";
            // cout << "residu gradb:  " << Linf(rhs_gradb[0]) / max_norm  << "  "<<endl;
            err_curr = max(Linf(rhs_gradb2[0]), max(Linf(rhs_g_ineq2[0]), max(Linf(rhs_g2[0]), max(Linf(rhs_rq2[0]), Linf(rhs_b2[0]))))) / max_norm;
            // cout << "residu:  " << err_curr << endl;
            if (i >= min_it_ref)
            {
                if (err_curr < 1e-6 || (error_prev > 0.0 && err_curr > 0.9 * error_prev))
                {
                    if (err_curr > 1e-6)
                    {
                        return -2;
                        // cout << "stopped it_ref because insufficient decrease err_curr:  " << err_curr << endl;
                    }
                    return 0;
                }
            }
            ls_->solve_rhs(
                &ocpkktmemory_,
                ux_test[0],
                lam_test[0],
                delta_s_test[0],
                sigma_total_cache[0],
                rhs_rq2[0],
                rhs_b2[0],
                rhs_g2[0],
                rhs_g_ineq2[0],
                rhs_gradb2[0]);
            // el = blasfeo_toc(&timer);
            // cout << "el time solveRHS " << el << endl;
            axpby(1.0, ux_test[0], 1.0, ux, ux);
            axpby(1.0, lam_test[0], 1.0, lam, lam);
            axpby(1.0, delta_s_test[0], 1.0, delta_s, delta_s);

            // prepare next iteration
            error_prev = err_curr;
        }
        printer_->level(1) << "WARNING: max number of refinement iterations reached, error: " << err_curr << endl;
    };
    return 0;
}
int FatropOCP::compute_scalings(
    double &obj_scale,
    FatropVecBF &x_scales,
    FatropVecBF &lam_scales,
    const FatropVecBF &grad_curr)
{
    return scaler_->compute_scalings(
        &ocpkktmemory_,
        obj_scale,
        x_scales,
        lam_scales,
        grad_curr);
};
int FatropOCP::eval_constraint_viol(
    const FatropVecBF &primal_vars,
    const FatropVecBF &slack_vars,
    FatropVecBF &constraint_violation)
{
    int res = ocp_->eval_contr_viol(
        &ocpkktmemory_,
        primal_vars,
        slack_vars,
        constraint_violation);
    return res;
};
int FatropOCP::eval_obj_grad(
    double obj_scale,
    const FatropVecBF &primal_vars,
    FatropVecBF &gradient)
{
    int res = ocp_->eval_obj_grad(
        &ocpkktmemory_,
        obj_scale,
        primal_vars,
        gradient);
    return res;
};
int FatropOCP::eval_obj(
    double obj_scale,
    const FatropVecBF &primal_vars,
    double &res)
{

    int resi = ocp_->eval_obj(
        &ocpkktmemory_,
        obj_scale,
        primal_vars,
        res);
    return resi;
};
int FatropOCP::eval_dual_inf(
    double obj_scale,
    const FatropVecBF &lam,
    const FatropVecBF &grad,
    FatropVecBF &du_inf)
{
    return duinfevaluator_.evaluate(
        &ocpkktmemory_,
        obj_scale,
        lam,
        grad,
        du_inf);
}
int FatropOCP::initialize_dual(
    const FatropVecBF &grad,
    FatropVecBF &dlam,
    // FatropVecBF &s_curr,
    const FatropVecBF &zL,
    const FatropVecBF &zU)
{
    // assume constraint jacobian evaluated
    OCPInitializer_.modify_kkt_ls_dual_estimate(&ocpkktmemory_, grad);
    s_dummy.SetConstant(0.0);
    s_zero.SetConstant(0.0);
    ux_dummy.SetConstant(0.0);
    sigma.SetConstant(1.0);
    axpy(-1.0, zL, zU, gradb);

    return ls_->solve_pd_sys(
        &ocpkktmemory_,
        0.0,
        0.0,
        ux_dummy,
        dlam,
        s_dummy,
        sigma,
        gradb);
    return 0;
}
int FatropOCP::initialize_slacks(FatropVecBF &s_curr)
{
    return OCPInitializer_.intialize_slack_variables(&ocpkktmemory_, s_curr);
}
NLPDims FatropOCP::get_nlp_dims() const
{
    return nlpdims_;
};
void FatropOCP::finalize()
{
}
void FatropOCP::reset()
{
    ocp_->reset();
}