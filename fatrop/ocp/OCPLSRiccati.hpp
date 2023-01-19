#ifndef OCPLSRICCATIINCLUDED
#define OCPLSRICCATIINCLUDED
#include "OCPKKT.hpp"
#include "OCPLinearSolver.hpp"
#include <cmath>
#define SUMMATION_ALG kahan_sum
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
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_s,
            const FatropVecBF &sigma_total,
            const FatropVecBF &gradb_total) override;
        // solve a KKT system
        int computeSDDeg(
            OCPKKTMemory *OCP,
            const double inertia_correction_w,
            const double inertia_correction_c,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_s,
            const FatropVecBF &sigma_total,
            const FatropVecBF &gradb_total);
        // solve a KKT system
        int
        computeSDnor(
            OCPKKTMemory *OCP,
            const double inertia_correction,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &delta_s,
            const FatropVecBF &sigma_total,
            const FatropVecBF &gradb_total);
        int GetRHS(
            OCPKKTMemory *OCP,
            const FatropVecBF &gradb_total,
            const FatropVecBF &rhs_rq,
            const FatropVecBF &rhs_b,
            const FatropVecBF &rhs_g,
            const FatropVecBF &rhs_g_ineq,
            const FatropVecBF &rhs_gradb);
        int ComputeMVProd(
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
            const FatropVecBF &rhs_gradb);
        int SolveRHS(
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
            const FatropVecBF &rhs_gradb);
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
        FatropMemoryVecBF rhs_rq;
        FatropMemoryVecBF rhs_b;
        FatropMemoryVecBF rhs_g;
        FatropMemoryVecBF rhs_g_ineq;
        FatropMemoryVecBF rhs_gradb;
        FatropMemoryVecBF rhs_rq2;
        FatropMemoryVecBF rhs_b2;
        FatropMemoryVecBF rhs_g2;
        FatropMemoryVecBF rhs_g_ineq2;
        FatropMemoryVecBF rhs_gradb2;
        FatropMemoryVecBF ux_test;
        FatropMemoryVecBF lam_test;
        FatropMemoryVecBF delta_s_test;
        FatropMemoryVecBF v_Ppt;
        FatropMemoryVecBF v_Hh;
        FatropMemoryVecBF v_AL;
        FatropMemoryVecBF v_RSQrqt_tilde;
        FatropMemoryVecBF v_Ggt_stripe;
        FatropMemoryVecBF v_Ggt_tilde;
        FatropMemoryVecBF v_GgLt;
        FatropMemoryVecBF v_RSQrqt_hat;
        FatropMemoryVecBF v_Llt;
        FatropMemoryVecBF v_Llt_shift;
        FatropMemoryVecBF v_GgIt_tilde;
        FatropMemoryVecBF v_GgLIt;
        FatropMemoryVecBF v_HhIt;
        FatropMemoryVecBF v_PpIt_hat;
        FatropMemoryVecBF v_LlIt;
        FatropMemoryVecBF v_Ggt_ineq_temp;
        MemoryPermMat Pl;
        MemoryPermMat Pr;
        MemoryPermMat PlI;
        MemoryPermMat PrI;
        FatropVector<int> gamma;
        FatropVector<int> rho;
        int rankI = 0;
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