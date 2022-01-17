#include "OCP/OCPTemplateBasic.hpp"
#include "OCP/OCPTemplate.hpp"
#include "OCP/OCPTemplAdapter.hpp"
#include "AUX/SmartPtr.hpp"
#include "AUX/FatropMemory.hpp"
#include "OCP/OCPAlg.hpp"
#include "OCP/OCPLinearSolver.hpp"
#include "OCP/OCPLSRiccati.hpp"
using namespace fatrop;
int main()
{
    MemoryAllocator fma;
    RefCountPtr<OCPTemplate> ocptemplatebasic =
        new OCPTemplateBasic(OCPTemplateBasic::from_shared_lib("./f.so", 10));
    RefCountPtr<OCP> ocptempladapter = new OCPTemplateAdapter(ocptemplatebasic, fma);
    RefCountPtr<OCPLinearSolver> ocplsriccati = new OCPLSRiccati(ocptempladapter->GetOCPDims(), fma);
    OCPAlg ocpalg(ocptempladapter, ocplsriccati, fma);
    fma.allocate();
}