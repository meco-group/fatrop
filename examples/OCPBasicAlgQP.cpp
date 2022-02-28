#include "ocp/BFOCPBasic.hpp"
// #include "ocp/BFOCP.hpp"
#include "ocp/BFOCPAdapter.hpp"
// #include "AUX/SmartPtr.hpp"
// #include "AUX/FatropMemory.hpp"
// #include "ocp/FatropOCP.hpp"
// #include "ocp/OCPLinearSolver.hpp"
#include "ocp/OCPLSRiccati.hpp"
// #include "AUX/FatropVector.hpp"
// #include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "ocp/OCPNoScaling.hpp"
#include "solver/FatropParams.hpp"
// #include "solver/FatropAlg.hpp"
// #ioclude "SPARSE/Sparseocp.hpp"
#include "solver/Filter.hpp"
#include "ocp/FatropOCP.hpp"
#include "solver/FatropAlg.hpp"
using namespace fatrop;
int main()
{
    // FatropMemoryMatBF test(2, 2, 1);
    // test[0].print();

    // RefCountPtr<BFOCP> ocptemplatebasic =
    //     new BFOCPBasic(BFOCPBasic::from_shared_lib("./../../OCP_specification/f.so", 100));
    // RefCountPtr<OCP> ocptempladapter = new BFOCPAdapter(ocptemplatebasic);
    // RefCountPtr<OCPLinearSolver> ocplsriccati = new OCPLSRiccati(ocptempladapter->GetOCPDims());
    // RefCountPtr<FatropParams> params = new FatropParams();
    // RefCountPtr<OCPScalingMethod> ocpscaler = new OCPNoScaling(params);
    // RefCountPtr<FatropNLP> fatropocp = new FatropOCP(ocptempladapter, ocplsriccati, ocpscaler);
    // RefCountPtr<FatropData> fatropdata = new FatropData(fatropocp->GetNLPDims(), params);
    // RefCountPtr<Filter> filter(new Filter(params->maxiter + 1));
    // RefCountPtr<Journaller> journaller(new Journaller(params->maxiter + 1));
    // RefCountPtr<LineSearch> linesearch = new BackTrackingLineSearch(params, fatropocp, fatropdata, filter, journaller);
    // RefCountPtr<FatropAlg> fatropalg = new FatropAlg(fatropocp, fatropdata, params, filter, linesearch, journaller);
    // blasfeo_timer timer;
    // VECSE(fatropdata->x_curr.nels(), 0.0, (VEC *)fatropdata->x_curr, 0);
    // VECSE(fatropdata->s_lower.nels(), 0.0, (VEC *)fatropdata->s_lower, 0);
    // VECSE(fatropdata->s_upper.nels(), INFINITY, (VEC *)fatropdata->s_upper, 0);
    // fatropdata->Initialize();
    // blasfeo_tic(&timer);
    // for (int i = 0; i < 1; i++)
    // {
    //     fatropalg->Optimize();
    // }
    // double el = blasfeo_toc(&timer);
    // cout << "el time " << el << endl;
    // // journaller->PrintIterations();
}