#ifndef NLPLSINCLUDED
#define NLPLSINCLUDED
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "aux/SmartPtr.hpp"
#include "OCPKKT.hpp"
namespace fatrop
{
    class OCPLinearSolver : public RefCountedObj
    {
    public:
        virtual int computeSD(
            OCPKKTMemory *OCP,
            const double intertia_correction_w,
            const double intertia_correction_c,
            const FatropVecBF &dprimal_vars,
            const FatropVecBF &dlam) = 0;
    };

} // namespace fatrop
#endif //  NLPLSINCLUDED