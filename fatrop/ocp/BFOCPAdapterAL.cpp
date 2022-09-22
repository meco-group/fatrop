#include "BFOCPAdapterAL.hpp"
using namespace fatrop;
int BFOCPAdapterAL::SetIneqsBounds(const vector<double> &lower_boundsin, const vector<double> &upper_boundsin)
{
    // copy(lower_boundsin, (ocptempl_->lower_bounds)[0]);
    // copy(upper_boundsin, (ocptempl_->upper_bounds)[0]);
    ocptempl_->lower_bounds = lower_boundsin;
    ocptempl_->upper_bounds = upper_boundsin;
    return 0;
}

int BFOCPAdapterAL::SetIneqLagrMult(const FatropVecBF &ineqlagrmultL, const FatropVecBF &ineqlagrmultU)
{
    copy(ineqlagrmultL, (ocptempl_->ineq_lagsL)[0]);
    copy(ineqlagrmultU, (ocptempl_->ineq_lagsU)[0]);
    return 0;
}

int BFOCPAdapterAL::SetPenalty(double penalty)
{
    ocptempl_->penalty = penalty;
    return 0;
}
int BFOCPAdapterAL::EvalInequalities(OCPKKTMemory *OCP,
                                     const FatropVecBF &primal_vars,
                                     FatropVecBF &g_ineq)
{
    const int K = ocptempl_->K;
    const int *offs_ux_p = (const int *)OCP->aux.ux_offs.data();
    int *offs_stageparams_p = (int *)offs_stageparams.data();
    double *stageparams_p = (double *)stageparams.data();
    double *globalparams_p = (double *)globalparams.data();
    int *nu_p = OCP->nu.data();
    int *offs_ineq = ocptempl_->ineqs_offsets.data();
    VEC *primal_vars_p = (VEC *)primal_vars;
    double *primal_data = primal_vars_p->pa;
    double *res_p = ((VEC *)g_ineq)->pa;
    for (int k = 0; k < K; k++)
    {
        const int offs = offs_ux_p[k];
        const int nu = nu_p[k];
        const double *inputs_k = primal_data + offs;
        const double *states_k = primal_data + (offs + nu);
        const double *stage_params_k = stageparams_p + offs_stageparams_p[k];
        const double *global_params_k = globalparams_p;
        double *res = res_p + offs_ineq[k];
        ocptempl_->eval_gineqk_AL(
            inputs_k,
            states_k,
            stage_params_k,
            global_params_k,
            res,
            k);
    }
    return 0;
}
int BFOCPAdapterAL::GetTotalNOIneqs()
{
    return ocptempl_->GetTotalNOIneqs();
}