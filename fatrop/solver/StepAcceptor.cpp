#include "solver/StepAcceptor.hpp"
using namespace fatrop;
StepAcceptor::StepAcceptor(
    const shared_ptr<FatropData> &fatropdata,
    const shared_ptr<Filter> &filter,
    const shared_ptr<FatropOptions> &fatropparams) : fatropdata_(fatropdata), filter_(filter), fatropoptions_(fatropparams)
{
    Initialize();
};
void StepAcceptor::Initialize()
{
}
void StepAcceptor::AcceptTrialStep()
{
}