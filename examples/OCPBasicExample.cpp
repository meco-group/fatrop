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
using namespace fatrop;
int main()
{
    RefCountPtr<BFOCP> ocptemplatebasic =
        new BFOCPBasic(BFOCPBasic::from_shared_lib("./f.so", 150));
    RefCountPtr<OCP> ocptempladapter = new BFOCPAdapter(ocptemplatebasic);
    RefCountPtr<OCPLinearSolver> ocplsriccati = new OCPLSRiccati(ocptempladapter->GetOCPDims());
    RefCountPtr<OCPScalingMethod> ocpscaler = new OCPNoScaling();
    FatropOCP ocpalg(ocptempladapter, ocplsriccati, ocpscaler);
    NLPDims nlpdims = ocpalg.GetNLPDims();
    blasfeo_timer timer;
    FatropMemoryVecBF ux(nlpdims.nvars, 2);
    FatropMemoryVecBF lags(nlpdims.neqs, 2);
    int N = 1000;
    VEC *prim_sc = (VEC *)ux[1];
    VEC *du_sc = (VEC *)lags[1];
    VECSE(nlpdims.nvars, 1.0, prim_sc, 0);
    VECSE(nlpdims.neqs, 1.0, du_sc, 0);
    blasfeo_tic(&timer);
    for (int i = 0; i < N; i++)
    {
        ocpalg.EvalHess(1.0, ux[0], ux[1], lags[0], lags[1]);
        ocpalg.EvalJac(ux[0], ux[1], lags[1]);
    }
    double el = blasfeo_toc(&timer);
    cout << "time elapsed " << el / N << endl;
}