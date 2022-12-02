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
            const double inertia_correction_w,
            const double inertia_correction_c,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_zL,
            const FatropVecBF &delta_zU,
            const FatropVecBF &delta_s,
            const FatropVecBF &sigma_L,
            const FatropVecBF &sigma_U,
            const FatropVecBF &gradb_L,
            const FatropVecBF &gradb_U,
            const FatropVecBF &gradb_plus,
            const FatropVecBF &zL_curr,
            const FatropVecBF &zU_curr) = 0;
    };

} // namespace fatrop
#endif //  NLPLSINCLUDED