#include "OCP/BFOCPBasic.hpp"
// #include "OCP/BFOCP.hpp"
#include "OCP/BFOCPAdapter.hpp"
// #include "AUX/SmartPtr.hpp"
// #include "AUX/FatropMemory.hpp"
// #include "OCP/FatropOCP.hpp"
// #include "OCP/OCPLinearSolver.hpp"
#include "OCP/OCPLSRiccati.hpp"
// #include "AUX/FatropVector.hpp"
// #include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "OCP/OCPNoScaling.hpp"
#include "SOLVER/FatropParams.hpp"
// #include "SOLVER/FatropAlg.hpp"
// #ioclude "SPARSE/SparseOCP.hpp"
#include "SOLVER/Filter.hpp"
#include "OCP/FatropOCP.hpp"
#include "SOLVER/FatropAlg.hpp"
using namespace fatrop;
int main()
{
    RefCountPtr<BFOCP> ocptemplatebasic =
        new BFOCPBasic(BFOCPBasic::from_shared_lib("./../../OCP_specification/f.so", 100));
    RefCountPtr<OCP> ocptempladapter = new BFOCPAdapter(ocptemplatebasic);
    RefCountPtr<OCPLinearSolver> ocplsriccati = new OCPLSRiccati(ocptempladapter->GetOCPDims());
    RefCountPtr<FatropParams> params = new FatropParams();
    RefCountPtr<OCPScalingMethod> ocpscaler = new OCPNoScaling(params);
    RefCountPtr<FatropNLP> fatropocp = new FatropOCP(ocptempladapter, ocplsriccati, ocpscaler);
    RefCountPtr<FatropData> fatropdata = new FatropData(fatropocp->GetNLPDims(), params);
    // VECSE(fatropdata->x_curr.nels(), 1.0, (VEC*)fatropdata->x_curr,0);
    RefCountPtr<Filter> filter(new Filter(params->maxiter + 1));
    RefCountPtr<LineSearch> linesearch = new BackTrackingLineSearch(params, fatropocp, fatropdata, filter);
    RefCountPtr<FatropAlg> fatropalg = new FatropAlg(fatropocp, fatropdata, params, filter,linesearch);
    blasfeo_timer timer;
    blasfeo_tic(&timer);

    fatropalg->Optimize();
    double el = blasfeo_toc(&timer);
    cout << "el time " << el << endl;
    }