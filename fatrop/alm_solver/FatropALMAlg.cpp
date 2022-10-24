#include "FatropALMAlg.hpp"
using namespace fatrop;
void FatropALMAlg::SetBounds(const vector<double> &lower_boundsin, const vector<double> &upper_boundsin)
{
    almdata_.lower_bounds = lower_boundsin;
    almdata_.upper_bounds = upper_boundsin;
};

void FatropALMAlg::SetInitial(const vector<double> &initial)
{
    almdata_.initial_x = initial;
};

void FatropALMAlg::GetSolution(vector<double> &sol)
{
    innersolver_.fatropdata_->x_curr.copyto(sol);
};

int FatropALMAlg::Optimize()
{
    int n_ineqs = almdata_.n_ineqs_;
    VECSE(n_ineqs,0.0, (VEC*) almdata_.auglags_Lcurr, 0);
    VECSE(n_ineqs,0.0, (VEC*) almdata_.auglags_Ucurr, 0);
    innersolver_.fatropdata_->x_initial = almdata_.initial_x;
    double penalty = 1e3;
    fatropnlpal_->SetPenalty(penalty);
    fatropnlpal_->SetIneqLagrMult(almdata_.auglags_Lcurr, almdata_.auglags_Ucurr);
    innersolver_.Optimize();
    for (int i = 0; i < 100; i++)
    {
        innersolver_.fatropdata_->x_initial.copy(innersolver_.fatropdata_->x_curr);
        fatropnlpal_->EvalInequalities(innersolver_.fatropdata_->x_curr, almdata_.ineq_curr);
        if(almdata_.MaxDist()<1e-3) return 0;
        almdata_.UpdateLags(penalty);
        fatropnlpal_->SetIneqLagrMult(almdata_.auglags_Lcurr, almdata_.auglags_Ucurr);
        if (true) penalty *= 1.2;
        fatropnlpal_->SetPenalty(penalty);
        innersolver_.Optimize();
    }
    return 0;
}
