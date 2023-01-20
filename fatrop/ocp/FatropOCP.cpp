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
                                                  gradb_total_cache(dims_.n_g_ineq_tot, 1),
                                                  sigma_total_cache(dims_.n_g_ineq_tot, 1){};
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

    return ls_->SolveRHS(
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
};
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