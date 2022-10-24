#ifndef OCPLSRICCATIINCLUDED
#define OCPLSRICCATIINCLUDED
#include "OCPKKT.hpp"
#include "OCPLinearSolver.hpp"
#include <cmath>
namespace fatrop
{
    bool check_reg(const int m, MAT *sA, const int ai, const int aj);
    class OCPLSRiccati : public OCPLinearSolver
    {
    public:
        OCPLSRiccati(const OCPDims &dims);
        // solve a KKT system
        int computeSD(
            OCPKKTMemory *OCP,
            const double inertia_correction_w,
            const double inertia_correction_c,
            const double mu,
            const double kappa_d,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &lam_curr,
            const FatropVecBF &s,
            const FatropVecBF &zL,
            const FatropVecBF &zU,
            const FatropVecBF &delta_zL,
            const FatropVecBF &delta_zU,
            const FatropVecBF &lower,
            const FatropVecBF &upper,
            const FatropVecBF &delta_s) override;
        // solve a KKT system
        int computeSDDeg(
            OCPKKTMemory *OCP,
            const double inertia_correction_w,
            const double inertia_correction_c,
            double mu,
            double kappa_d,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &lam_curr,
            const FatropVecBF &s,
            const FatropVecBF &zL,
            const FatropVecBF &zU,
            const FatropVecBF &delta_zL,
            const FatropVecBF &delta_zU,
            const FatropVecBF &lower,
            const FatropVecBF &upper,
            const FatropVecBF &delta_s);
        // solve a KKT system
        int
        SolveInitialization(
            OCPKKTMemory *OCP,
            const FatropVecBF &lam,
            const FatropVecBF &ux_dummy,
            const FatropVecBF &s_dummy,
            const FatropVecBF &zL,
            const FatropVecBF &zU,
            const FatropVecBF &lower,
            const FatropVecBF &upper) override;
        // solve a KKT system
        int
        computeSDnor(
            OCPKKTMemory *OCP,
            const double inertia_correction,
            const double mu,
            const double kappa_d,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &lam_curr,
            const FatropVecBF &s,
            const FatropVecBF &zL,
            const FatropVecBF &zU,
            const FatropVecBF &delta_zL,
            const FatropVecBF &delta_zU,
            const FatropVecBF &lower,
            const FatropVecBF &upper,
            const FatropVecBF &delta_s);
        int computeResidual(
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
            const FatropVecBF &delta_s,
            FatropVecBF &residual) ;
        FatropMemoryMatBF Ppt;
        FatropMemoryMatBF Hh;
        FatropMemoryMatBF AL;
        FatropMemoryMatBF RSQrqt_tilde;
        FatropMemoryMatBF Ggt_stripe;
        FatropMemoryMatBF Ggt_tilde;
        FatropMemoryMatBF GgLt;
        FatropMemoryMatBF RSQrqt_hat;
        FatropMemoryMatBF Llt;
        FatropMemoryMatBF Llt_shift; // needed because feature not implemented yet
        FatropMemoryMatBF GgIt_tilde;
        FatropMemoryMatBF GgLIt;
        FatropMemoryMatBF HhIt;
        FatropMemoryMatBF PpIt_hat;
        FatropMemoryMatBF LlIt;
        FatropMemoryMatBF Ggt_ineq_temp;
        MemoryPermMat Pl;
        MemoryPermMat Pr;
        MemoryPermMat PlI;
        MemoryPermMat PrI;
        FatropVector<int> gamma;
        FatropVector<int> rho;
        struct LastUsed
        {
            int rankI = 0;
            double inertia_correction = 0;
            double kappa_d = 0;
            double mu = 0;
        } lastused_;
    };
};     // namespace
#endif // OCPLSRICCATIINCLUDED