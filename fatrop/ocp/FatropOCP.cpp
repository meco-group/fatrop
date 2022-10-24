#include "ocp/FatropOCP.hpp"
using namespace fatrop;
FatropOCP::FatropOCP(
    const shared_ptr<OCP> &ocp,
    const shared_ptr<OCPLinearSolver> &ls,
    const shared_ptr<OCPScalingMethod> &scaler) : ocp_(ocp), ls_(ls), scaler_(scaler), ocpkktmemory_(ocp_->GetOCPDims()){};
int FatropOCP::EvalHess(
    double obj_scale,
    const FatropVecBF &primal_vars,
    const FatropVecBF &lam) 
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    int res = ocp_->evalHess(
        &ocpkktmemory_,
        obj_scale,
        primal_vars,
        lam);
    hess_time += blasfeo_toc(&timer);
    hess_count += 1;
    return res;
};
int FatropOCP::EvalJac(
    const FatropVecBF &primal_vars,
    const FatropVecBF &slack_vars) 
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    int res = ocp_->evalJac(
        &ocpkktmemory_,
        primal_vars,
        slack_vars);
    jac_time += blasfeo_toc(&timer);
    jac_count += 1;
    return res;
};
int FatropOCP::ComputeSD(
    const double inertia_correction_w,
    const double inertia_correction_c,
    const double mu,
    const double kappa_d,
    const FatropVecBF &dprimal_vars,
    const FatropVecBF &dlam,
    const FatropVecBF &lam_curr,
    const FatropVecBF &s,
    const FatropVecBF &zL_curr,
    const FatropVecBF &zU_curr,
    const FatropVecBF &delta_zL,
    const FatropVecBF &delta_zU,
    const FatropVecBF &lower_bound,
    const FatropVecBF &upper_bound,
    const FatropVecBF &delta_s) 
{
    // ls_ = RefCountPtr<OCPLinearSolver>(new Sparse_OCP(ocp_->GetOCPDims(), ocpkktmemory_));
    return ls_->computeSD(
        &ocpkktmemory_,
        inertia_correction_w,
        inertia_correction_c,
        mu,
        kappa_d,
        dprimal_vars,
        dlam,
        lam_curr,
        s,
        zL_curr,
        zU_curr,
        delta_zL,
        delta_zU,
        lower_bound,
        upper_bound,
        delta_s);
};
int FatropOCP::ComputeResidual(
    const double inertia_correction_w,
    const double inertia_correction_c,
    const double mu,
    const double kappa_d,
    const FatropVecBF &dprimal_vars,
    const FatropVecBF &dlam,
    const FatropVecBF &lam_curr,
    const FatropVecBF &s,
    const FatropVecBF &zL_curr,
    const FatropVecBF &zU_curr,
    const FatropVecBF &delta_zL,
    const FatropVecBF &delta_zU,
    const FatropVecBF &lower_bound,
    const FatropVecBF &upper_bound,
    const FatropVecBF &delta_s, 
    FatropVecBF &residual) 
{
    // ls_ = RefCountPtr<OCPLinearSolver>(new Sparse_OCP(ocp_->GetOCPDims(), ocpkktmemory_));
    return ls_->computeResidual(
        &ocpkktmemory_,
        inertia_correction_w,
        inertia_correction_c,
        mu,
        kappa_d,
        dprimal_vars,
        dlam,
        lam_curr,
        s,
        zL_curr,
        zU_curr,
        delta_zL,
        delta_zU,
        lower_bound,
        upper_bound,
        delta_s, 
        residual);
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
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    int res = ocp_->EvalConstraintViolation(
        &ocpkktmemory_,
        primal_vars,
        slack_vars,
        constraint_violation);
    cv_time += blasfeo_toc(&timer);
    cv_count += 1;
    return res;
};
int FatropOCP::EvalGrad(
    double obj_scale,
    const FatropVecBF &primal_vars,
    FatropVecBF &gradient) 
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    int res = ocp_->EvalGrad(
        &ocpkktmemory_,
        obj_scale,
        primal_vars,
        gradient);
    grad_time += blasfeo_toc(&timer);
    grad_count += 1;
    return res;
};
int FatropOCP::EvalObj(
    double obj_scale,
    const FatropVecBF &primal_vars,
    double &res) 
{

    blasfeo_timer timer;
    blasfeo_tic(&timer);
    int resi = ocp_->EvalObj(
        &ocpkktmemory_,
        obj_scale,
        primal_vars,
        res);
    obj_time += blasfeo_toc(&timer);
    obj_count += 1;
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
    const FatropVecBF &ux_dummy,
    const FatropVecBF &s_dummy,
    FatropVecBF &s_curr,
    const FatropVecBF &zL,
    const FatropVecBF &zU,
    const FatropVecBF &lower,
    const FatropVecBF &upper) 
{
    // assume constraint jacobian evaluated
    OCPInitializer_.AdaptKKTInitial(&ocpkktmemory_, grad, s_curr);
    return ls_->SolveInitialization(&ocpkktmemory_, dlam, ux_dummy, s_dummy, zL, zU, lower, upper);
}
NLPDims FatropOCP::GetNLPDims() const 
{
    NLPDims res;
    // states + inputs
    res.nvars = sum(ocpkktmemory_.nu) + sum(ocpkktmemory_.nx);
    // stagewise equality contraints + dynamics constraints + slack constraints
    res.neqs = sum(ocpkktmemory_.ng) + sum(ocpkktmemory_.nx) - ocpkktmemory_.nx.at(0) + sum(ocpkktmemory_.ng_ineq);
    res.nineqs = sum(ocpkktmemory_.ng_ineq);
    return res;
};
void FatropOCP::Finalize() 
{
    FE_time = hess_time + jac_time + cv_time + grad_time + obj_time;
    cout << "hess time: " << hess_time << ", count: " << hess_count << endl;
    cout << "jac time:  " << jac_time << ", count: " << jac_count << endl;
    cout << "cv time:   " << cv_time << ", count: " << cv_count << endl;
    cout << "grad time: " << grad_time << ", count: " << grad_count << endl;
    cout << "obj time:  " << obj_time << ", count: " << obj_count << endl;
    cout << "FE time:   " << FE_time << endl;
}
void FatropOCP::Reset() 
{
    FE_time = 0.0;
    hess_time = 0.0;
    jac_time = 0.0;
    cv_time = 0.0;
    grad_time = 0.0;
    obj_time = 0.0;
    hess_count = 0;
    jac_count = 0;
    cv_count = 0;
    grad_count = 0;
    obj_count = 0;
}