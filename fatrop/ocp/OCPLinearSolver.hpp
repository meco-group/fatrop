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
            const FatropVecBF &delta_s,
            const FatropVecBF &sigma_total,
            const FatropVecBF &gradb_total) = 0;
        virtual int ComputeMVProd(
            OCPKKTMemory *OCP,
            const double inertia_correction_w,
            const double inertia_correction_c,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_s,
            const FatropVecBF &sigma_total,
            const FatropVecBF &rhs_rq,
            const FatropVecBF &rhs_b,
            const FatropVecBF &rhs_g,
            const FatropVecBF &rhs_g_ineq,
            const FatropVecBF &rhs_gradb_total) = 0;
        virtual int SolveRHS(
            OCPKKTMemory *OCP,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_s,
            const FatropVecBF &sigma_total,
            const FatropVecBF &rhs_rq,
            const FatropVecBF &rhs_b,
            const FatropVecBF &rhs_g,
            const FatropVecBF &rhs_g_ineq,
            const FatropVecBF &rhs_gradb_total) = 0;
        virtual int GetRHS(
            OCPKKTMemory *OCP,
            const FatropVecBF &gradb_total,
            const FatropVecBF &rhs_rq,
            const FatropVecBF &rhs_b,
            const FatropVecBF &rhs_g,
            const FatropVecBF &rhs_g_ineq,
            const FatropVecBF &rhs_gradb) = 0;
    };

} // namespace fatrop
#endif //  NLPLSINCLUDED