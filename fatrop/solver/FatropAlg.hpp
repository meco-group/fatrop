#ifndef FATROPALGINCLUDED
#define FATROPALGINCLUDED
#include "templates/NLPAlg.hpp"
#include "FatropData.hpp"
#include "Filter.hpp"
#include "LineSearch.hpp"
#include <cmath>
#include "IterationData.hpp"
#include <memory>
#include "FatropStats.hpp"
#include <limits>
#include <solver/FatropPrinter.hpp>
// #include "AlgorithmQuantities.hpp"
#ifdef ENABLE_MULTITHREADING
#include "aux/Worker.hpp"
#endif

namespace fatrop
{
    // TODO: get rid of FatropApplication
    class FatropAlg 
    {
    public:
        FatropAlg(
            const std::shared_ptr<FatropNLP> &fatropnlp,
            const std::shared_ptr<FatropData> &fatropdata,
            const std::shared_ptr<FatropOptions> &fatropparams,
            const std::shared_ptr<Filter> &filter,
            const std::shared_ptr<LineSearch> &linesearch,
            const std::shared_ptr<Journaller> &journaller);
        void initialize() ;
        void reset() ;
        void set_bounds(const std::vector<double> &lower, const std::vector<double> &upper) ;
        void set_initial(const std::vector<double> &initial) ;
        void get_solution(std::vector<double> &sol) ;
        int optimize() ;
        int eval_lag_hess();
        int eval_constr_viol();
        inline int eval_constr_viol_curr();
        inline int eval_constr_viol_trial();
        inline int eval_obj_grad_curr();
        double eval_objective_curr();
        double eval_objective_trial();
        int eval_dual_infeasiblity();
        inline int perform_initializiation();
        int solve_pd_sys(double inertia_correction_w, double inertia_correction_c, double mu);
        std::shared_ptr<FatropNLP> fatropnlp_;
        std::shared_ptr<FatropData> fatropdata_;
        std::shared_ptr<FatropOptions> fatropoptions_;
        std::shared_ptr<Filter> filter_;
        std::shared_ptr<LineSearch> linesearch_;
        std::shared_ptr<Journaller> journaller_;
        FatropStats get_stats() 
        {
            return stats;
        };

    public:
        double tol;
        double acceptable_tol;
        int acceptable_iter;
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
        double theta_min;
        int max_watchdog_steps;
        bool warm_start_init_point;
        // bool first_try_watchdog;
        FatropStats stats;
    };
} // namespace fatrop
#endif // FATROPALGINCLUDED