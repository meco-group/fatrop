#ifndef BFOCPADAPTERALINCLUDED
#define BFOCPADAPTERALINCLUDED
#include "BFOCPAdapter.hpp"
#include "BFOCPAL.hpp"
#include "OCPAL.hpp"
// but contains methods to set parameters specifically to inner problem (bounds, lagrangemultipliers, ...)
namespace fatrop
{
    class BFOCPAdapterAL : public BFOCPAdapter, OCPAL
    {
    public:
        BFOCPAdapterAL(const shared_ptr<BFOCPAL> &ocptempl) : BFOCPAdapter(ocptempl), ocptempl_(ocptempl)
        {
        }
        int SetIneqsBounds(const FatropVecBF &lower_boundsin, const FatropVecBF &upper_boundsin);
        int SetIneqLagrMult(const FatropVecBF &ineqlagrmult);
        int SetPenalty(double penalty);
        const shared_ptr<BFOCPAL> ocptempl_;
    };
} // namespace fatrop
#endif // BFOCPADAPTERALINCLUDED