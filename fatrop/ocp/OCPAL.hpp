#ifndef OCPALINCLUDED
#define OCPALINCLUDED
#include "OCP.hpp"
namespace fatrop
{
    class OCPAL : public OCP
    {
    public:
        virtual int SetIneqsBounds(const FatropVecBF &lower_boundsin, const FatropVecBF &upper_boundsin) = 0;
        virtual int SetIneqLagrMult(const FatropVecBF &ineqlagrmult) = 0;
        virtual int SetPenalty(double penalty) = 0;
    };
} // namespace fatrop
#endif // OCPALINCLUDED