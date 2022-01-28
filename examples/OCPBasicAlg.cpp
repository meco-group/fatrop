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
// #include "SPARSE/SparseOCP.hpp"
#include "OCP/FatropOCP.hpp"
#include "SOLVER/FatropAlg.hpp"
using namespace fatrop;
int main()
{
    RefCountPtr<BFOCP> ocptemplatebasic =
        new BFOCPBasic(BFOCPBasic::from_shared_lib("./../../OCP_specification/f.so", 3));
    RefCountPtr<OCP> ocptempladapter = new BFOCPAdapter(ocptemplatebasic);
    RefCountPtr<OCPLinearSolver> ocplsriccati = new OCPLSRiccati(ocptempladapter->GetOCPDims());
    RefCountPtr<FatropParams> params = new FatropParams();
    RefCountPtr<OCPScalingMethod> ocpscaler = new OCPNoScaling(params);
    RefCountPtr<FatropNLP> fatropocp = new FatropOCP(ocptempladapter, ocplsriccati, ocpscaler);
    RefCountPtr<FatropData> fatropdata = new FatropData(fatropocp->GetNLPDims(), params);
    RefCountPtr<FatropAlg> fatropalg = new FatropAlg(fatropocp, fatropdata, params);
    fatropalg->Optimize();
}