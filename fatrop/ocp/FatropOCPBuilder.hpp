#ifndef FATROPOCPBUILDERINCLUDED
#define FATROPOCPBUILDERINCLUDED
#include "ocp/OCP.hpp"
#include "ocp/FatropOCP.hpp"
#include "ocp/OCPAbstract.hpp"
#include "ocp/OCPAdapter.hpp"
#include "ocp/OCPLSRiccati.hpp"
#include "ocp/OCPNoScaling.hpp"
namespace fatrop
{
    class FatropOCPBuilder
    {
    public:
        FatropOCPBuilder(const shared_ptr<OCPAbstract> &ocp, const shared_ptr<FatropOptions> &fatropparams) : ocp_(ocp), fatropoptions_(fatropparams)
        {
        }
        shared_ptr<FatropOCP> Build()
        {
            shared_ptr<OCPAdapter> adapter = make_shared<OCPAdapter>(ocp_);
            return Build(adapter);
        }

        shared_ptr<FatropOCP> Build(shared_ptr<OCPAdapter> &adapter)
        {
            return make_shared<FatropOCP>(adapter, make_shared<OCPLSRiccati>(adapter->GetOCPDims(), fatropoptions_), make_shared<OCPNoScaling>(fatropoptions_), fatropoptions_);
        }

    private:
        const shared_ptr<OCPAbstract> ocp_;
        const shared_ptr<FatropOptions> fatropoptions_;
    };
}

#endif // !OCPALGBUILDERINCLUDED