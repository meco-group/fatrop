#ifndef OCPEVALUATORINCLUDED
#define OCPEVALUATORINCLUDED
#include "OCPKKT.hpp"
#include "BFOCP.hpp"
#include "solver/FatropData.hpp"
#include "aux/SmartPtr.hpp"
#include "OCP.hpp"
#include <memory>
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP->name)
#define AUXMACRO(type, name, suffix) type name##suffix = ((type)OCP->aux.name)
#define SOLVERMACRO(type, name, suffix) type name##suffix = ((type)name)
// #ifdef ENABLE_MULTITHREADING
// #include <omp.h>
// #endif

// Example: you can make a for loop parallel with the following code. But to be further investigated to do this in an efficient way...
// #ifdef ENABLE_MULTITHREADING
// #pragma omp parallel for
// <your for loop>
// #endif

namespace fatrop
{
    using namespace std;
    class BFOCPAdapter : public OCP // public OCP -> also include KKTmemory, OCPDims, ...
    {
    public:
        BFOCPAdapter(const shared_ptr<BFOCP> &ocptempl_) : K(ocptempl_->get_horizon_length()),
                                                           nuexpr(TransformRange<int>(0, K, [&ocptempl_](int k)
                                                                                      { return ocptempl_->get_nuk(k); })),
                                                           nxexpr(TransformRange<int>(0, K, [&ocptempl_](int k)
                                                                                      { return ocptempl_->get_nxk(k); })),
                                                           ngexpr(TransformRange<int>(0, K, [&ocptempl_](int k)
                                                                                      { return ocptempl_->get_ngk(k); })),
                                                           ngineqexpr(TransformRange<int>(0, K, [&ocptempl_](int k)
                                                                                          { return ocptempl_->get_ng_ineq_k(k); })),
                                                           nstageparamsexpr(TransformRange<int>(0, K, [&ocptempl_](int k)
                                                                                          { return ocptempl_->get_n_stage_params_k(k); })), offs_stageparams(offsets(nstageparamsexpr)), stageparams(sum(nstageparamsexpr), 0.0), globalparams(ocptempl_->get_n_global_parmas(), 0.0), ocptempl(ocptempl_)
        {
            x_dummy = vector<double>(max(nxexpr), 0.0);
        }
        int evalHess(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &lam) override;
        int evalJac(
            OCPKKTMemory *OCP,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars) override;
        int EvalConstraintViolation(
            OCPKKTMemory *OCP,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            FatropVecBF &constraint_violation) override;
        int EvalGrad(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            FatropVecBF &gradient) override;
        int EvalObj(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            double &res);
        int EvalDynamics(
            OCPKKTMemory *OCP,
            const int k,
            const FatropVecBF &uk,
            const FatropVecBF &xk,
            FatropVecBF &xkp1);
        OCPDims GetOCPDims() const override
        {
            return OCPDims(ocptempl->get_horizon_length(), nuexpr, nxexpr, ngexpr, ngineqexpr);
        }
    public:
        void SetParams(const vector<double> &stage_params_in, const vector<double> &global_params_in) override;
        void SetInitial(const shared_ptr<FatropData> &fatropdata, vector<double> &initial_u, vector<double> &initial_x) override;
        void GetSolution(const shared_ptr<FatropData> &fatropdata, vector<double> &u, vector<double> &x) override;
        double& GetGlobalParams()
        {
            return globalparams.at(0);
        }
        double& GetStageParams()
        {
            return stageparams.at(0);
        }

    public:
        int K;
        FatropVector<int> nuexpr;
        FatropVector<int> nxexpr;
        FatropVector<int> ngexpr;
        FatropVector<int> ngineqexpr;
        FatropVector<int> nstageparamsexpr;
        FatropVector<int> offs_stageparams;
        vector<double> stageparams;
        vector<double> globalparams;
        vector<double> x_dummy;
    private:
        shared_ptr<BFOCP> ocptempl;
    };
} // namespace fatrop

#endif // OCPEVALUATORINCLUDED