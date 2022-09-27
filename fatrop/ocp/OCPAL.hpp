#ifndef OCPALINCLUDED
#define OCPALINCLUDED
#include "OCP.hpp"
namespace fatrop
{
    class OCPAL : public OCP
    {
    public:
        virtual int SetIneqsBounds(const vector<double> &lower_boundsin, const vector<double> &upper_boundsin) = 0;
        virtual int SetIneqLagrMult(const FatropVecBF &ineqlagrmultL, const FatropVecBF &ineqlagrmultU) = 0;
        virtual int ResetIneqLagrMult() = 0;
        virtual int SetPenalty(double penalty) = 0;
        virtual int EvalInequalities(OCPKKTMemory *OCP,
                                     const FatropVecBF &primal_vars,
                                     FatropVecBF &g_ineq) = 0;
        virtual int GetTotalNOIneqs() = 0;
    };
} // namespace fatrop
#endif // OCPALINCLUDED