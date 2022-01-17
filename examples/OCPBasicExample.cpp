#include "OCP/OCPTemplateBasic.hpp"
#include "OCP/OCPTemplate.hpp"
#include "OCP/OCPTemplAdapter.hpp"
#include "AUX/SmartPtr.hpp"
#include "AUX/FatropMemory.hpp"
#include "OCP/OCPAlg.hpp"
using namespace fatrop;
int main()
{
    MemoryAllocator fma;
    RefCountPtr<OCPTemplate> ocptemplatebasic =
        new OCPTemplateBasic(OCPTemplateBasic::from_shared_lib("./f.so", 10));
    OCPTemplateAdapter ocptempladapter(ocptemplatebasic, fma);
}