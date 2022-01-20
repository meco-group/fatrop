#ifndef OCPALGINCLUDED
#define OCPALGINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "TEMPLATES/NLPAlg.hpp"
#include "AUX/SmartPtr.hpp"
#include "OCPKKT.hpp"
#include "OCPLinearSolver.hpp"
#include "AUX/FatropMemory.hpp"
namespace fatrop
{
    class FatropOCP : public FatropNLP
    {
    public:
        FatropOCP(const RefCountPtr<OCP> &ocp, const RefCountPtr<OCPLinearSolver> &ls) : ocp_(ocp), ls_(ls), ocpkktmemory_(ocp_->GetOCPDims()){};
        int EvalHess(
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &scales_primal_vars,
            const FatropVecBF &lam,
            const FatropVecBF &scales_lam)
        {
            return ocp_->evalHess(
                &ocpkktmemory_,
                obj_scale,
                primal_vars,
                scales_primal_vars,
                lam,
                scales_lam);
        };
        int EvalJac(
            const FatropVecBF &primal_vars,
            const FatropVecBF &scales_primal_vars,
            const FatropVecBF &scales_lam)
        {
            return ocp_->evalJac(
                &ocpkktmemory_,
                primal_vars,
                scales_primal_vars,
                scales_lam);
        };
        int ComputeSD(
            const double inertia_correction,
            const FatropVecBF &dprimal_vars,
            const FatropVecBF &dlam)
        {
            return ls_->computeSD(
                &ocpkktmemory_,
                inertia_correction,
                dprimal_vars,
                dlam);
        };

    public:
        RefCountPtr<OCP> ocp_;
        RefCountPtr<OCPLinearSolver> ls_;
        OCPKKTMemory ocpkktmemory_;
    };
} // namespace fatrop
#endif //  OCPALGINCLUDED