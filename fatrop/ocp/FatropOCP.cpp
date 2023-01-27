#include "ocp/FatropOCP.hpp"
using namespace fatrop;
FatropOCP::FatropOCP(
    const shared_ptr<OCP> &ocp,
    const shared_ptr<OCPLinearSolver> &ls,
    const shared_ptr<OCPScalingMethod> &scaler) : ocp_(ocp), dims_(ocp_->GetOCPDims()),
                                                  nlpdims_({sum(dims_.nx + dims_.nu), sum(dims_.ng + dims_.ng_ineq + dims_.nx) - dims_.nx.at(0), sum(dims_.ng_ineq)}), ls_(ls), scaler_(scaler), ocpkktmemory_(dims_), s_memvec(nlpdims_.nineqs, 4), ux_memvec(nlpdims_.nvars, 1),
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
}
int FatropOCP::EvalHess(
    double obj_scale,
    const FatropVecBF &primal_vars,
    const FatropVecBF &lam)
{
    int res = ocp_->evalHess(
        &ocpkktmemory_,
        obj_scale,
        primal_vars,
        lam);
    return res;
};
int FatropOCP::EvalJac(
    const FatropVecBF &primal_vars,
    const FatropVecBF &slack_vars)
{
    int res = ocp_->evalJac(
        &ocpkktmemory_,
        primal_vars,
        slack_vars);
    return res;
};
int FatropOCP::ComputeSD(
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

    return ls_->computeSD(
        &ocpkktmemory_,
        inertia_correction_w,
        inertia_correction_c,
        ux,
        lam,
        delta_s,
        sigma_total,
        gradb_total);
};

int FatropOCP::SolveSOC(
    const FatropVecBF &ux,
    const FatropVecBF &lam,
    const FatropVecBF &delta_s,
    const FatropVecBF &constraint_violation)
{
    bool it_ref = true;
    /// todo avoid retrieving unnecessary rhs'es
    ls_->GetRHS(
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

    ls_->SolveRHS(
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
        double max_norm = std::max(Linf(rhs_gradb[0]), std::max(Linf(rhs_g_ineq[0]), std::max(Linf(rhs_g[0]), std::max(Linf(rhs_rq[0]), Linf(rhs_b[0])))));
        double error_prev = -1.0;
        for (int i = 0; i < 5; i++)
        {
            ls_->ComputeMVProd(
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
            err_curr = std::max(Linf(rhs_gradb2[0]), std::max(Linf(rhs_g_ineq2[0]), std::max(Linf(rhs_g2[0]), std::max(Linf(rhs_rq2[0]), Linf(rhs_b2[0]))))) / max_norm;
            // cout << "residu:  " << err_curr << endl;
            if (err_curr < 1e-8 || (error_prev > 0.0 && err_curr > 0.9 * error_prev))
            {
                if (err_curr > 1e-12)
                {
                    // cout << "stopped it_ref because insufficient decrease err_curr:  " << err_curr << endl;
                }
                return 0;
            }
            ls_->SolveRHS(
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
        cout << "WARNING: max number of refinement iterations reached, error: " << err_curr << endl;
    };
    return 0;
}
int FatropOCP::ComputeScalings(
    double &obj_scale,
    FatropVecBF &x_scales,
    FatropVecBF &lam_scales,
    const FatropVecBF &grad_curr)
{
    return scaler_->ComputeScalings(
        &ocpkktmemory_,
        obj_scale,
        x_scales,
        lam_scales,
        grad_curr);
};
int FatropOCP::EvalConstraintViolation(
    const FatropVecBF &primal_vars,
    const FatropVecBF &slack_vars,
    FatropVecBF &constraint_violation)
{
    int res = ocp_->EvalConstraintViolation(
        &ocpkktmemory_,
        primal_vars,
        slack_vars,
        constraint_violation);
    return res;
};
int FatropOCP::EvalGrad(
    double obj_scale,
    const FatropVecBF &primal_vars,
    FatropVecBF &gradient)
{
    int res = ocp_->EvalGrad(
        &ocpkktmemory_,
        obj_scale,
        primal_vars,
        gradient);
    return res;
};
int FatropOCP::EvalObj(
    double obj_scale,
    const FatropVecBF &primal_vars,
    double &res)
{

    int resi = ocp_->EvalObj(
        &ocpkktmemory_,
        obj_scale,
        primal_vars,
        res);
    return resi;
};
int FatropOCP::EvalDuInf(
    double obj_scale,
    const FatropVecBF &lam,
    const FatropVecBF &grad,
    FatropVecBF &du_inf)
{
    return duinfevaluator_.DuInfEval(
        &ocpkktmemory_,
        obj_scale,
        lam,
        grad,
        du_inf);
}
int FatropOCP::Initialization(
    const FatropVecBF &grad,
    FatropVecBF &dlam,
    FatropVecBF &s_curr,
    const FatropVecBF &zL,
    const FatropVecBF &zU)
{
    // assume constraint jacobian evaluated
    OCPInitializer_.AdaptKKTInitial(&ocpkktmemory_, grad, s_curr);
    s_dummy.SetConstant(0.0);
    s_zero.SetConstant(0.0);
    ux_dummy.SetConstant(0.0);
    sigma.SetConstant(1.0);
    axpy(-1.0, zL, zU, gradb);

    return ls_->computeSD(
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
NLPDims FatropOCP::GetNLPDims() const
{
    return nlpdims_;
};
void FatropOCP::Finalize()
{
}
void FatropOCP::Reset()
{
}