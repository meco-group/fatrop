#include "OCP/BFOCPBasic.hpp"
#include "OCP/BFOCP.hpp"
#include "OCP/BFOCPAdapter.hpp"
#include "AUX/SmartPtr.hpp"
#include "AUX/FatropMemory.hpp"
#include "OCP/OCPAlg.hpp"
#include "OCP/OCPLinearSolver.hpp"
#include "OCP/OCPLSRiccati.hpp"
using namespace fatrop;
int main()
{
    MemoryAllocator fma;
    RefCountPtr<BFOCP> ocptemplatebasic =
        new BFOCPBasic(BFOCPBasic::from_shared_lib("./f.so", 10));
    RefCountPtr<OCP> ocptempladapter = new BFOCPAdapter(ocptemplatebasic, fma);
    RefCountPtr<OCPLinearSolver> ocplsriccati = new OCPLSRiccati(ocptempladapter->GetOCPDims(), fma);
    FatropOCP ocpalg(ocptempladapter, ocplsriccati, fma);
    fma.allocate();
}