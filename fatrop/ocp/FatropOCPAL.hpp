#ifndef FATROPOCPALINCLUDED
#define FATROPOCPALINCLUDED
#include "FatropOCP.hpp"
#include "templates/FatropNLPAL.hpp"
#include "ocp/OCPAL.hpp"
using namespace std;
namespace fatrop
{
    class FatropOCPAL : public FatropOCP, FatropNLPAL
    {
        FatropOCPAL(
            const shared_ptr<OCPAL> &ocp,
            const shared_ptr<OCPLinearSolver> &ls,
            const shared_ptr<OCPScalingMethod> &scaler) : FatropOCP(ocp, ls, scaler), ocp_(ocp)
        {
        }
            int SetIneqsBounds(const FatropVecBF &lower_boundsin, const FatropVecBF &upper_boundsin)
        {
            return ocp_->SetIneqsBounds(lower_boundsin, upper_boundsin);
        };
        int SetIneqLagrMult(const FatropVecBF &ineqlagrmult)
        {
            return ocp_->SetIneqLagrMult(ineqlagrmult);
        };
        int SetPenalty(double penalty)
        {
            return ocp_->SetPenalty(penalty);
        };
        const shared_ptr<OCPAL> ocp_;
    };

    }  // namespace fatrop
#endif // FATROPOCPALINCLUDED