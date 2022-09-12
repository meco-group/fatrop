#ifndef LINESEARCHDDPINCLUDED
#define LINESEARCHDDPINCLUDED
#include "solver/LineSearch.hpp"
#include "ocp/OCPLSRiccati.hpp"
#include "ocp/OCP.hpp"
namespace fatrop
{
    // the LineSearchDDP is just an adapation of the BacktrackingLinesearch
    class LineSearchDDP : public BackTrackingLineSearch
    {
    public:
        LineSearchDDP(
            const shared_ptr<FatropParams> &fatropparams,
            const shared_ptr<FatropNLP> &nlp,
            const shared_ptr<FatropData> &fatropdata,
            const shared_ptr<Filter> &filter,
            const shared_ptr<Journaller> &journaller,
            const shared_ptr<OCPLSRiccati>& ocplsriccati,
            OCPKKTMemory* OCPmem,
            const shared_ptr<OCP>& ocpinterface)
            : BackTrackingLineSearch(fatropparams, nlp, fatropdata, filter, journaller), ocplsriccati_(ocplsriccati), OCP_(OCPmem), ocpinterface_(ocpinterface){};
        int TryStep(double alpha_pr, double alpha_du) const override;
    private:
        const shared_ptr<OCPLSRiccati> ocplsriccati_;
        OCPKKTMemory* OCP_;
        const shared_ptr<OCP> ocpinterface_;
    };

} // namespace fatrop
#endif // LINESEARCHDDPINCLUDED