#ifndef FATROPOCPBUILDERINCLUDED
#define FATROPOCPBUILDERINCLUDED
#include "ocp/OCP.hpp"
#include "ocp/FatropOCP.hpp"
#include "ocp/OCPAbstract.hpp"
#include "ocp/BFOCPAdapter.hpp"
#include "ocp/OCPLSRiccati.hpp"
#include "ocp/OCPNoScaling.hpp"
namespace fatrop
{
    class FatropOCPBuilder
    {
    public:
        FatropOCPBuilder(const shared_ptr<OCPAbstract> &ocp, const shared_ptr<FatropParams>& fatropparams ) : ocp_(ocp), fatropparams_(fatropparams)
        {
        }
        shared_ptr<FatropOCP> Build()
        {
           shared_ptr<BFOCPAdapter> adapter(make_shared<BFOCPAdapter>(ocp_));
           return make_shared<FatropOCP>(adapter, make_shared<OCPLSRiccati>(adapter->GetOCPDims()), make_shared<OCPNoScaling>(fatropparams_)); 
        }

    private:
        const shared_ptr<OCPAbstract> ocp_;
        const shared_ptr<FatropParams> fatropparams_;
    };
}

#endif // !OCPALGBUILDERINCLUDED