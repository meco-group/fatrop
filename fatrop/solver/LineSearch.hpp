#ifndef LINESEARCHINCLUDED
#define LINESEARCHINCLUDED
#include "AlgStrategy.hpp"
#include "IterationData.hpp"
#include "templates/NLPAlg.hpp"
#include "solver/FatropData.hpp"
#include "solver/Filter.hpp"
#include "solver/FatropPrinter.hpp"
#include <cmath>
#include <memory>
using namespace std;
namespace fatrop
{
    struct LineSearchInfo
    {
        int ls = 0;
        bool first_rejected_by_filter = false;
        bool last_rejected_by_filter = false;
    };
    class LineSearch : public AlgStrategy
    {
    public:
        LineSearch(
            const shared_ptr<FatropParams> &fatropparams,
            const shared_ptr<FatropNLP> &nlp,
            const shared_ptr<FatropData> &fatropdata);
        virtual LineSearchInfo FindAcceptableTrialPoint(double mu, bool small_sd, bool from_backup) = 0;
        inline int EvalCVNext();
        double EvalObjNext();
        void Reset();
        virtual int TryStep(double alpha_pr, double alpha_du) const;
        virtual int InitSoc() const;
        virtual int ExitSoc() const;
        virtual int CalcSoc(double alpha) const;
        shared_ptr<FatropNLP> fatropnlp_;
        shared_ptr<FatropData> fatropdata_;
        int eval_cv_count;
        int eval_obj_count;
        int eval_cv_time;
        int eval_obj_time;
    };

    class BackTrackingLineSearch : public LineSearch
    {
    public:
        BackTrackingLineSearch(
            const shared_ptr<FatropParams> &fatropparams,
            const shared_ptr<FatropNLP> &nlp,
            const shared_ptr<FatropData> &fatropdata,
            const shared_ptr<Filter> &filter,
            const shared_ptr<Journaller> &journaller);
        void Initialize();
        LineSearchInfo FindAcceptableTrialPoint(double mu, bool small_sd, bool from_backup);
        shared_ptr<Filter> filter_;
        shared_ptr<Journaller> journaller_;
        double s_phi;
        double delta;
        double s_theta;
        double gamma_theta;
        double gamma_phi;
        double eta_phi;
        double gamma_alpha;
    };
} // namespace fatrop
#endif // LINESEARCHINCLUDED