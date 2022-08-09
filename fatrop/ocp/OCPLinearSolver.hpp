#ifndef NLPLSINCLUDED
#define NLPLSINCLUDED
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "aux/SmartPtr.hpp"
#include "OCPKKT.hpp"
namespace fatrop
{
    class OCPLinearSolver 
    {
    public:
        virtual int computeSD(
            OCPKKTMemory *OCP,
            const double intertia_correction_w,
            const double intertia_correction_c,
            const double mu,
            const double kappa_d,
            const FatropVecBF &dprimal_vars,
            const FatropVecBF &dlam,
            const FatropVecBF &lam_curr,
            const FatropVecBF &s,
            const FatropVecBF &zL,
            const FatropVecBF &zU,
            const FatropVecBF &delta_zL,
            const FatropVecBF &delta_zU,
            const FatropVecBF &lower,
            const FatropVecBF &upper,
            const FatropVecBF &delta_s) = 0;

        virtual int
        SolveInitialization(
            OCPKKTMemory *OCP,
            const FatropVecBF &lam,
            const FatropVecBF &ux_dummy,
            const FatropVecBF &s_dummy,
            const FatropVecBF &zL,
            const FatropVecBF &zU,
            const FatropVecBF &lower,
            const FatropVecBF &upper) = 0;
    };

} // namespace fatrop
#endif //  NLPLSINCLUDED