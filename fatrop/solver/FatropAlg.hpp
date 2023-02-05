#ifndef FATROPALGINCLUDED
#define FATROPALGINCLUDED
#include "aux/SmartPtr.hpp"
#include "templates/NLPAlg.hpp"
#include "FatropData.hpp"
#include "Filter.hpp"
#include "LineSearch.hpp"
#include "StepAcceptor.hpp"
#include <cmath>
#include "IterationData.hpp"
#include "templates/FatropApplication.hpp"
#include <memory>
#include "FatropStats.hpp"
#include <limits>
using namespace std;
// #include "AlgorithmQuantities.hpp"
#ifdef ENABLE_MULTITHREADING
#include "aux/Worker.hpp"
#endif

namespace fatrop
{
    class FatropAlg : public FatropApplication
    {
    public:
        FatropAlg(
            const shared_ptr<FatropNLP> &fatropnlp,
            const shared_ptr<FatropData> &fatropdata,
            const shared_ptr<FatropParams> &fatropparams,
            const shared_ptr<Filter> &filter,
            const shared_ptr<LineSearch> &linesearch,
            const shared_ptr<Journaller> &journaller);
        void Initialize() override;
        void Reset() override;
        void SetBounds(const vector<double> &lower, const vector<double> &upper) override;
        void SetInitial(const vector<double> &initial) override;
        void GetSolution(vector<double> &sol) override;
        int Optimize() override;
        inline int EvalHess();
        inline int EvalJac();
        inline int EvalCVCurr();
        inline int EvalCVNext();
        inline int EvalGradCurr();
        double EvalObjCurr();
        double EvalObjNext();
        int EvalDuInf();
        inline int Initialization();
        int ComputeSD(double inertia_correction_w, double inertia_correction_c, double mu);
        void WarmStart() override
        {
            fatropdata_->x_initial.copy(fatropdata_->x_curr);
        };
        shared_ptr<FatropNLP> fatropnlp_;
        shared_ptr<FatropData> fatropdata_;
        shared_ptr<FatropParams> fatropparams_;
        shared_ptr<Filter> filter_;
        shared_ptr<LineSearch> linesearch_;
        shared_ptr<Journaller> journaller_;
        FatropStats GetStats() 
        {
            return stats;
        };

    public:
        double tol;
        double acceptable_tol;
        double acceptable_iter;
        int maxiter;

    private:
        double lammax;
        double mu0;
        double kappa_eta;
        double kappa_mu;
        double theta_mu;
        double delta_w0;
        double delta_wmin;
        double kappa_wmin;
        double kappa_wplus;
        double kappa_wplusem;
        double delta_c_stripe;
        double kappa_c;
        double kappa_d;
        int max_watchdog_steps;
        // bool first_try_watchdog;
        FatropStats stats;
    };
} // namespace fatrop
#endif // FATROPALGINCLUDED