#ifndef OCPINCLUDED
#define OCPINCLUDED
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "OCPDims.hpp"
#include "aux/SmartPtr.hpp"
#include "ocp/OCPKKT.hpp"
#include "solver/FatropData.hpp"
namespace fatrop
{ 
    /** \brief interface class for OCP operations*/
    class OCP : public RefCountedObj
    {
    public:
        virtual int evalHess(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &lam) = 0;
        virtual int evalJac(
            OCPKKTMemory *OCP,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars) = 0;
        virtual int EvalConstraintViolation(
            OCPKKTMemory *OCP,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            FatropVecBF &constraint_violation) = 0;
        virtual int EvalGrad(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            FatropVecBF &gradient) = 0;
        virtual int EvalObj(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            double &res) = 0;
        virtual OCPDims GetOCPDims() const = 0;
        virtual void SetParams(const vector<double> &stage_params_in, const vector<double> &global_params_in) = 0;
        virtual void SetInitial(const int K, const RefCountPtr<FatropData> &fatropdata, vector<double> &initial_u, vector<double> &initial_x) = 0;
    };
} // namespace fatrop
#endif // OCPINCLUDED