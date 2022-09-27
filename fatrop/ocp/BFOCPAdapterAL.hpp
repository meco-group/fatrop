#ifndef BFOCPADAPTERALINCLUDED
#define BFOCPADAPTERALINCLUDED
#include "BFOCPAdapter.hpp"
#include "BFOCPAL.hpp"
#include "OCPAL.hpp"
// but contains methods to set parameters specifically to inner problem (bounds, lagrangemultipliers, ...)
namespace fatrop
{
    class BFOCPAdapterAL : public OCPAL, public BFOCPAdapter
    {
    public:
        BFOCPAdapterAL(const shared_ptr<BFOCPAL> &ocptempl) : BFOCPAdapter(ocptempl), ocptempl_(ocptempl)
        {
        }
        int evalHess(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &lam)
        {
            return BFOCPAdapter::evalHess(
                OCP,
                obj_scale,
                primal_vars,
                lam);
        };
        int evalJac(
            OCPKKTMemory *OCP,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars)
        {
            return BFOCPAdapter::evalJac(
                OCP,
                primal_vars,
                slack_vars);
        };
        int EvalConstraintViolation(
            OCPKKTMemory *OCP,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            FatropVecBF &constraint_violation)
        {
            return BFOCPAdapter::EvalConstraintViolation(
                OCP,
                primal_vars,
                slack_vars,
                constraint_violation);
        };
        int EvalGrad(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            FatropVecBF &gradient)
        {
            return BFOCPAdapter::EvalGrad(
                OCP,
                obj_scale,
                primal_vars,
                gradient);
        };
        int EvalObj(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            double &res)
        {
            return BFOCPAdapter::EvalObj(
                OCP,
                obj_scale,
                primal_vars,
                res);
        };
        int EvalDynamics(
            OCPKKTMemory *OCP,
            const int k,
            const FatropVecBF &uk,
            const FatropVecBF &xk,
            FatropVecBF &xkp1)
        {
            return BFOCPAdapter::EvalDynamics(
                OCP,
                k,
                uk,
                xk,
                xkp1);
        };
        OCPDims GetOCPDims() const override
        {
            return BFOCPAdapter::GetOCPDims();
        }
        void SetParams(const vector<double> &stage_params_in, const vector<double> &global_params_in) override
        {
            BFOCPAdapter::SetParams(stage_params_in, global_params_in);
        };
        void SetInitial(const int K, const shared_ptr<FatropData> &fatropdata, vector<double> &initial_u, vector<double> &initial_x)
        {
            BFOCPAdapter::SetInitial(K, fatropdata, initial_u, initial_x);
        };
        int SetIneqsBounds(const vector<double> &lower_boundsin, const vector<double> &upper_boundsin) override;
        int SetIneqLagrMult(const FatropVecBF &ineqlagrmultL, const FatropVecBF &ineqlagrmultU) override;
        int ResetIneqLagrMult() override;
        int SetPenalty(double penalty) override;
        int EvalInequalities(OCPKKTMemory *OCP,
                             const FatropVecBF &primal_vars,
                             FatropVecBF &g_ineq) override;
        int GetTotalNOIneqs();
        const shared_ptr<BFOCPAL> ocptempl_;
    };
} // namespace fatrop
#endif // BFOCPADAPTERALINCLUDED