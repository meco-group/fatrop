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
        FatropOCPBuilder(const std::shared_ptr<OCPAbstract> &ocp, const std::shared_ptr<FatropOptions> &fatropparams) : ocp_(ocp), fatropoptions_(fatropparams)
        {
        }
        std::shared_ptr<FatropOCP> build()
        {
            std::shared_ptr<OCPAdapter> adapter = std::make_shared<OCPAdapter>(ocp_);
            return build(adapter);
        }

        std::shared_ptr<FatropOCP> build(std::shared_ptr<OCPAdapter> &adapter)
        {
            return std::make_shared<FatropOCP>(adapter, std::make_shared<OCPLSRiccati>(adapter->get_ocp_dims(), fatropoptions_), std::make_shared<OCPNoScaling>(fatropoptions_), fatropoptions_);
        }

    private:
        const std::shared_ptr<OCPAbstract> ocp_;
        const std::shared_ptr<FatropOptions> fatropoptions_;
    };
}

#endif // !OCPALGBUILDERINCLUDED