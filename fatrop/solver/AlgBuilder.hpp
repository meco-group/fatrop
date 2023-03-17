#ifndef ALBBUILDERINCLUDED
#define ALBBUILDERINCLUDED
// #include "NLPAlg.hpp"
#include "FatropAlg.hpp"
namespace fatrop
{
    class AlgBuilder
    {
    public:
        void BuildFatropAlgObjects(const shared_ptr<FatropNLP> &nlp,
                                    const shared_ptr<FatropOptions> &fatropparams,
                                   shared_ptr<FatropData> &fatropdata,
                                   shared_ptr<Journaller> &journaller)
        {
            fatropdata = make_shared<FatropData>(nlp->GetNLPDims(), fatropparams);
            journaller = make_shared<Journaller>(fatropparams->maxiter + 1);
            fatropdata_ = fatropdata; // keep this around for building the algorithm
            journaller_ = journaller;
            nlp_ = nlp;
            fatropoptions_ = fatropparams;
        }
        shared_ptr<FatropAlg> BuildAlgorithm()
        {
            // TODO unsafe if maxiter is changed during application
            shared_ptr<Filter> filter = make_shared<Filter>(fatropoptions_->maxiter + 1);
            shared_ptr<LineSearch> linesearch = make_shared<BackTrackingLineSearch>(fatropoptions_, nlp_, fatropdata_, filter, journaller_);
            return make_shared<FatropAlg>(nlp_, fatropdata_, fatropoptions_, filter, linesearch, journaller_);
        }
    private:
        shared_ptr<FatropNLP> nlp_;
        shared_ptr<FatropOptions> fatropoptions_;
        shared_ptr<FatropData> fatropdata_;
        shared_ptr<Journaller> journaller_;
    };
};

#endif // ALBBUILDERINCLUDED