#ifndef ALBBUILDERINCLUDED
#define ALBBUILDERINCLUDED
// #include "NLPAlg.hpp"
#include "FatropAlg.hpp"
namespace fatrop
{
    class AlgBuilder
    {
    public:
        void BuildFatropAlgObjects(const std::shared_ptr<FatropNLP> &nlp,
                                    const std::shared_ptr<FatropOptions> &fatropparams,
                                   std::shared_ptr<FatropData> &fatropdata,
                                   std::shared_ptr<Journaller> &journaller)
        {
            fatropdata = std::make_shared<FatropData>(nlp->GetNLPDims(), fatropparams);
            journaller = std::make_shared<Journaller>(fatropparams->maxiter + 1);
            fatropdata_ = fatropdata; // keep this around for building the algorithm
            journaller_ = journaller;
            nlp_ = nlp;
            fatropoptions_ = fatropparams;
        }
        std::shared_ptr<FatropAlg> BuildAlgorithm()
        {
            // TODO unsafe if maxiter is changed during application
            std::shared_ptr<Filter> filter = std::make_shared<Filter>(fatropoptions_->maxiter + 1);
            std::shared_ptr<LineSearch> linesearch = std::make_shared<BackTrackingLineSearch>(fatropoptions_, nlp_, fatropdata_, filter, journaller_);
            return std::make_shared<FatropAlg>(nlp_, fatropdata_, fatropoptions_, filter, linesearch, journaller_);
        }
    private:
        std::shared_ptr<FatropNLP> nlp_;
        std::shared_ptr<FatropOptions> fatropoptions_;
        std::shared_ptr<FatropData> fatropdata_;
        std::shared_ptr<Journaller> journaller_;
    };
};

#endif // ALBBUILDERINCLUDED