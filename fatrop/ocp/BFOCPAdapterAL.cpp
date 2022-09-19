#include "BFOCPAdapterAL.hpp"
using namespace fatrop;
int BFOCPAdapterAL::SetIneqsBounds(const FatropVecBF &lower_boundsin, const FatropVecBF &upper_boundsin)
{
    copy(lower_boundsin, (ocptempl_->lower_bounds)[0]);
    copy(upper_boundsin, (ocptempl_->upper_bounds)[0]);
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