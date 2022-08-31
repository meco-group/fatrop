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
        BFOCPAdapter(const shared_ptr<BFOCP> &ocptempl_) : nuexpr(shared_ptr<BFOCP>(ocptempl_)), nxexpr(shared_ptr<BFOCP>(ocptempl_)), ngexpr(shared_ptr<BFOCP>(ocptempl_)), ngineqexpr(shared_ptr<BFOCP>(ocptempl_)), nstageparamsexpr(shared_ptr<BFOCP>(ocptempl_)), offs_stageparams(offsets(nstageparamsexpr)), stageparams(sum(nstageparamsexpr), 0.0), globalparams(ocptempl_->get_n_global_parmas(), 0.0), ocptempl(ocptempl_)
        {
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
        OCPDims GetOCPDims() const override
        {
            return OCPDims(ocptempl->get_horizon_length(), nuexpr, nxexpr, ngexpr, ngineqexpr);
        }

    private:
        class nxExpr : public VecExpr<nxExpr, int>
        {
        public:
            nxExpr(const shared_ptr<BFOCP> &parent) : parent(parent){};
            int getEl(const int ai) const { return parent->get_nxk(ai); };
            int size() const { return parent->get_horizon_length(); };

        private:
            const shared_ptr<BFOCP> parent;
        };
        class nuExpr : public VecExpr<nuExpr, int>
        {
        public:
            nuExpr(const shared_ptr<BFOCP> &parent) : parent(parent){};
            int getEl(const int ai) const { return parent->get_nuk(ai); };
            int size() const { return parent->get_horizon_length(); };

        private:
            const shared_ptr<BFOCP> parent;
        };
        class ngExpr : public VecExpr<ngExpr, int>
        {
        public:
            ngExpr(const shared_ptr<BFOCP> &parent) : parent(parent){};
            int getEl(const int ai) const { return parent->get_ngk(ai); };
            int size() const { return parent->get_horizon_length(); };

        private:
            const shared_ptr<BFOCP> parent;
        };
        class ngIneqExpr : public VecExpr<ngIneqExpr, int>
        {
        public:
            ngIneqExpr(const shared_ptr<BFOCP> &parent) : parent(parent){};
            int getEl(const int ai) const { return parent->get_ng_ineq_k(ai); };
            int size() const { return parent->get_horizon_length(); };

        private:
            const shared_ptr<BFOCP> parent;
        };
        class nStageParamsExpr : public VecExpr<nStageParamsExpr, int>
        {
        public:
            nStageParamsExpr(const shared_ptr<BFOCP> &parent) : parent(parent){};
            int getEl(const int ai) const { return parent->get_n_stage_params_k(ai); };
            int size() const { return parent->get_horizon_length(); };

        private:
            const shared_ptr<BFOCP> parent;
        };

    public:
        void SetParams(const vector<double> &stage_params_in, const vector<double> &global_params_in) override;
        void SetInitial(const int K, const shared_ptr<FatropData> &fatropdata, vector<double> &initial_u, vector<double> &initial_x);

    public:
        nuExpr nuexpr;
        nxExpr nxexpr;
        ngExpr ngexpr;
        ngIneqExpr ngineqexpr;
        nStageParamsExpr nstageparamsexpr;
        FatropVector<int> offs_stageparams;
        vector<double> stageparams;
        vector<double> globalparams;

    private:
        shared_ptr<BFOCP> ocptempl;
    };
} // namespace fatrop

#endif // OCPEVALUATORINCLUDED