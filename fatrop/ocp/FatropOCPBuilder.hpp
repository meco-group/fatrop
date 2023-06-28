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
        FatropOCPBuilder(const std::shared_ptr<OCPAbstract> &ocp, const std::shared_ptr<FatropOptions> &fatropparams, const std::shared_ptr<FatropPrinter> &printer) : ocp_(ocp), fatropoptions_(fatropparams), fatropprinter_(printer)
        {
        }
        std::shared_ptr<FatropOCP> build()
        {
            std::shared_ptr<OCPAdapter> adapter = std::make_shared<OCPAdapter>(ocp_);
            return build(adapter);
        }

        std::shared_ptr<FatropOCP> build(std::shared_ptr<OCPAdapter> &adapter)
        {
            return std::make_shared<FatropOCP>(adapter, std::make_shared<OCPLSRiccati>(adapter->get_ocp_dims(), fatropoptions_, fatropprinter_), std::make_shared<OCPNoScaling>(fatropoptions_), fatropoptions_, fatropprinter_);
        }

    private:
        const std::shared_ptr<OCPAbstract> ocp_;
        const std::shared_ptr<FatropOptions> fatropoptions_;
        const std::shared_ptr<FatropPrinter> fatropprinter_;
    };
}

#endif // !OCPALGBUILDERINCLUDED