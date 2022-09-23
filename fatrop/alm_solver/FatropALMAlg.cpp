#include "FatropALMAlg.hpp"
using namespace fatrop;
void FatropALMAlg::SetBounds(const vector<double> &lower_boundsin, const vector<double> &upper_boundsin) 
        {
            almdata_.lower_bounds = lower_boundsin;
            almdata_.upper_bounds = upper_boundsin;
        };

void FatropALMAlg::SetInitial(const vector<double>& initial) 
        {
            innersolver_.fatropdata_->x_curr = initial;
        };

void FatropALMAlg::GetSolution(vector<double>& sol) 
        {
            innersolver_.fatropdata_->x_curr.copyto(sol);
        };

int FatropALMAlg::Optimize() 
        {
            double penalty = 1e0;
            fatropnlpal_->SetPenalty(penalty);
            innersolver_.Optimize();
            for (int i = 0; i < 20; i++)
            {
                innersolver_.fatropdata_->x_initial.copy(innersolver_.fatropdata_->x_curr);
                fatropnlpal_->EvalInequalities(innersolver_.fatropdata_->x_curr, almdata_.ineq_curr);
                almdata_.UpdateLags(penalty);
                fatropnlpal_->SetIneqLagrMult(almdata_.auglags_Lcurr, almdata_.auglags_Ucurr);
                fatropnlpal_->SetPenalty(penalty);
                innersolver_.Optimize();
                if(i<5)
                penalty *= 10;
            }
            return 0;
        }

