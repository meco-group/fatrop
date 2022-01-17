#ifndef OCPALGINCLUDED
#define OCPALGINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "TEMPLATES/NLPAlg.hpp"
#include "AUX/SmartPtr.hpp"
#include "OCPKKT.hpp"
#include "AUX/FatropMemory.hpp"
namespace fatrop
{
    class OCPAlg : public NLPAlg
    {
    public:
        OCPAlg(const RefCountPtr<OCP> &ocp, MemoryAllocator& fma) : ocp_(ocp), ocpkktmemory_(ocp_->GetOCPDims(), fma){};

    private:
        OCPKKTMemory ocpkktmemory_;
        RefCountPtr<OCP> ocp_;
    };
} // namespace fatrop
#endif //  OCPALGINCLUDED