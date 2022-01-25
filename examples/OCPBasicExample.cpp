#include "OCP/BFOCPBasic.hpp"
#include "OCP/BFOCP.hpp"
#include "OCP/BFOCPAdapter.hpp"
#include "AUX/SmartPtr.hpp"
#include "AUX/FatropMemory.hpp"
#include "OCP/FatropOCP.hpp"
#include "OCP/OCPLinearSolver.hpp"
#include "OCP/OCPLSRiccati.hpp"
#include "AUX/FatropVector.hpp"
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "OCP/OCPNoScaling.hpp"
#include "SOLVER/FatropParams.hpp"
#include "SOLVER/FatropAlg.hpp"

using namespace fatrop;
int main()
{
    RefCountPtr<BFOCP> ocptemplatebasic =
        new BFOCPBasic(BFOCPBasic::from_shared_lib("./f.so", 150));
    RefCountPtr<OCP> ocptempladapter = new BFOCPAdapter(ocptemplatebasic);
    RefCountPtr<OCPLinearSolver> ocplsriccati = new OCPLSRiccati(ocptempladapter->GetOCPDims());
    RefCountPtr<FatropParams> params = new FatropParams();
    RefCountPtr<OCPScalingMethod> ocpscaler = new OCPNoScaling(params);
    FatropOCP ocpalg(ocptempladapter, ocplsriccati, ocpscaler);
    NLPDims nlpdims = ocpalg.GetNLPDims();
    blasfeo_timer timer;
    FatropMemoryVecBF ux(nlpdims.nvars, 4);
    FatropMemoryVecBF lags(nlpdims.neqs, 4);
    int N = 1000;
    FatropVecBF scalesx = ux[1];
    FatropVecBF dux = ux[3];
    FatropVecBF grad = ux[2];
    FatropVecBF cv = lags[2];
    FatropVecBF scaleslam = lags[1];
    FatropVecBF dlam = lags[3];
    // FatropAlg fatropalg(ocpalg)
    double obj_sc = 1.0;
    ocpalg.ComputeScalings(
        obj_sc,
        scalesx,
        scaleslam,
        scaleslam);
    blasfeo_tic(&timer);
    for (int i = 0; i < N; i++)
    {
        ocpalg.EvalConstraintViolation(ux[0], ux[1], scaleslam, cv);
        ocpalg.EvalGrad(1.0, ux[0], ux[1], grad);
        ocpalg.EvalHess(1.0, ux[0], ux[1], lags[0], lags[1]);
        ocpalg.EvalJac(ux[0], ux[1], lags[1]);
        ocpalg.ComputeSD(0.0, dux, dlam);
    }
    double el = blasfeo_toc(&timer);
    cout << "time elapsed " << el / N << endl;
}