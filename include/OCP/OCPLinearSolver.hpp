#ifndef NLPLSINCLUDED
#define NLPLSINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "AUX/SmartPtr.hpp"
#include "OCPKKT.hpp"
namespace fatrop
{
    class OCPLinearSolver : public RefCountedObj
    {
    public:
        virtual int computeSD(
            OCPKKTMemory *OCP,
            const double intertia_correction,
            const FatropVecBF &dprimal_vars,
            const FatropVecBF &dlam) = 0;
    };

} // namespace fatrop
#endif //  NLPLSINCLUDED