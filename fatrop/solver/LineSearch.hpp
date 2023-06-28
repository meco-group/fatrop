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
            const std::shared_ptr<FatropOptions> &fatropparams,
            const std::shared_ptr<FatropNLP> &nlp,
            const std::shared_ptr<FatropData> &fatropdata, const std::shared_ptr<FatropPrinter> &printer);
        virtual LineSearchInfo find_acceptable_trial_point(double mu, bool small_sd, bool from_backup) = 0;
        inline int eval_constr_viol_trial();
        double eval_obj_trial();
        void reset();
        virtual int update_trial_step(double alpha_pr, double alpha_du) const;
        virtual int initialize_second_order_correction() const;
        virtual int exit_second_order_correction() const;
        virtual int compute_second_order_correction(double alpha) const;
        std::shared_ptr<FatropNLP> fatropnlp_;
        std::shared_ptr<FatropData> fatropdata_;
        std::shared_ptr<FatropPrinter> printer_;
        int eval_cv_count;
        int eval_obj_count;
        double eval_cv_time;
        double eval_obj_time;
    };

    class BackTrackingLineSearch : public LineSearch
    {
    public:
        BackTrackingLineSearch(
            const std::shared_ptr<FatropOptions> &fatropparams,
            const std::shared_ptr<FatropNLP> &nlp,
            const std::shared_ptr<FatropData> &fatropdata,
            const std::shared_ptr<Filter> &filter,
            const std::shared_ptr<Journaller> &journaller, const std::shared_ptr<FatropPrinter> &printer);
        void initialize();
        LineSearchInfo find_acceptable_trial_point(double mu, bool small_sd, bool from_backup);
        std::shared_ptr<Filter> filter_;
        std::shared_ptr<Journaller> journaller_;
        double s_phi;
        double delta;
        double s_theta;
        double gamma_theta;
        double gamma_phi;
        double eta_phi;
        double gamma_alpha;
        bool accept_every_trial_step = false;
        int max_soc;
    };
} // namespace fatrop
#endif // LINESEARCHINCLUDED