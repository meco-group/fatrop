#ifndef FATROPALGINCLUDED
#define FATROPALGINCLUDED
#include "AUX/SmartPtr.hpp"
#include "TEMPLATES/NLPAlg.hpp"
#include "FatropData.hpp"
#include "Filter.hpp"
#include "LineSearch.hpp"
#include "StepAcceptor.hpp"
#include <cmath>
#include "IterationData.hpp"
// #include "AlgorithmQuantities.hpp"
namespace fatrop
{
    class FatropAlg : public RefCountedObj
    {
    public:
        FatropAlg(
            const RefCountPtr<FatropNLP> &fatropnlp,
            const RefCountPtr<FatropData> &fatropdata,
            const RefCountPtr<FatropParams> &fatropparams,
            const RefCountPtr<Filter> &filter,
            const RefCountPtr<LineSearch> &linesearch,
            const RefCountPtr<Journaller> &journaller)
            : fatropnlp_(fatropnlp),
              fatropdata_(fatropdata),
              fatropparams_(fatropparams),
              filter_(filter),
              linesearch_(linesearch),
              journaller_(journaller)
        {
            Initialize();
        }
        void Initialize()
        {
            lammax = fatropparams_->lammax;
            maxiter = fatropparams_->maxiter;
            tol = fatropparams_->tol;
            mu0 = fatropparams_->mu0;
            kappa_eta = fatropparams_->kappa_eta;
            kappa_mu = fatropparams_->kappa_mu;
            theta_mu = fatropparams_->theta_mu;
            delta_w0 = fatropparams_->delta_w0;
            delta_wmin = fatropparams_->delta_wmin;
            kappa_wmin = fatropparams_->kappa_wmin;
            kappa_wplus = fatropparams_->kappa_wplus;
            kappa_wplusem = fatropparams_->kappa_wplusem;
            delta_c_stripe = fatropparams_->delta_c_stripe;
            kappa_c = fatropparams_->kappa_c;
            // todo avoid reallocation when maxiter doesn't change
            // filter_ = RefCountPtr<Filter>(new Filter(maxiter + 1));
        }
        void Reset()
        {
            filter_->Reset();
            fatropdata_->Reset();
            journaller_->Reset();
        }
        int Optimize()
        {
            Reset();
            // optimization algorithm parameters
            const double mu_min = tol / 10;
            // optimization variables
            double mu = mu0;
            double delta_w_last = 0.0;
            // initialization
            EvalJac(); // todo twice evaluation
            EvalGradCurr();
            Initialization();
            if (fatropdata_->LamLinfCalc() < lammax)
            {
                // cout << "accepted lam " << endl;
                fatropdata_->AcceptInitialization();
            }
            EvalCVCurr();
            fatropdata_->theta_min = fatropparams_->theta_min * MAX(1.0, fatropdata_->CVL1Curr());
            int ls = 0;
            double deltaw = 0;
            double deltac = 0.0;
            for (int i = 0; i < maxiter; i++)
            {
                fatropdata_->obj_curr = EvalObjCurr();
                EvalJac();      // needed for dual inf
                EvalGradCurr(); // needed for dual inf
                EvalDuInf();
                IterationData &it_curr = journaller_->it_curr;
                it_curr.iter = i;
                it_curr.mu = mu;
                it_curr.objective = EvalObjCurr();
                it_curr.constraint_violation = fatropdata_->CVLinfCurr();
                it_curr.du_inf = fatropdata_->DuInfLinfCurr();
                it_curr.ls = ls;
                it_curr.reg = deltaw;
                journaller_->Push();
                journaller_->PrintIterations();
                double emu = fatropdata_->EMuCurr(0.0);
                if (emu < tol)
                {
                    cout << "found solution :) " << endl;
                    return 0;
                }
                // update mu
                // todo make a seperate class
                while (mu > mu_min && (fatropdata_->EMuCurr(mu) <= kappa_eta * mu))
                {
                    mu = MAX(mu_min, MIN(kappa_mu * mu, pow(mu, theta_mu)));
                    filter_->Reset();
                }
                // Hessian is necessary for calculating search direction
                EvalHess();
                // todo make an update class for regularization
                double deltac_candidate = delta_c_stripe * pow(mu, kappa_c);
                deltaw = 0.0;
                deltac = 0.0;
                int regularity = ComputeSD(deltaw, deltac);
                if (regularity < 0)
                {
                    deltac = deltac_candidate;
                    regularity = ComputeSD(deltaw, deltac);
                    cout << "Jac degenerate" << endl;
                }
                int increase_counter = 0;
                if (regularity > 0) // regularization is necessary
                {
                    deltaw = (delta_w_last == 0.0) ? delta_w0 : MAX(delta_wmin, kappa_wmin * delta_w_last);
                    regularity = ComputeSD(deltaw, deltac);
                    if ((deltac == 0.0) && (regularity < 0))
                    {
                        deltac = deltac_candidate;
                        regularity = ComputeSD(deltaw, deltac);
                        cout << "Jac degenerate" << endl;
                    }
                    while (regularity > 0)
                    {
                        increase_counter++;
                        deltaw = (delta_w_last == 0.0) ? kappa_wplusem * deltaw : kappa_wplus * deltaw;
                        regularity = ComputeSD(deltaw, deltac);
                        if ((deltac == 0.0) && (regularity < 0))
                        {
                            deltac = deltac_candidate;
                            regularity = ComputeSD(deltaw, deltac);
                            cout << "Jac degenerate" << endl;
                        }
                    }
                    delta_w_last = deltaw;
                }
                // cout << "regularization  " << (deltaw) << endl;
                // cout << "step size " << Linf(fatropdata_->delta_x) << endl;
                ls = linesearch_->FindAcceptableTrialPoint();
            }
            return 0;
        }
        inline int EvalHess()
        {
            return fatropnlp_->EvalHess(
                fatropdata_->obj_scale,
                fatropdata_->x_curr,
                fatropdata_->lam_curr);
        }
        inline int EvalJac()
        {
            return fatropnlp_->EvalJac(
                fatropdata_->x_curr);
        }
        inline int EvalCVCurr()
        {
            return fatropnlp_->EvalConstraintViolation(
                fatropdata_->x_curr,
                fatropdata_->g_curr);
        }
        inline int EvalGradCurr()
        {
            return fatropnlp_->EvalGrad(
                fatropdata_->obj_scale,
                fatropdata_->x_curr,
                fatropdata_->grad_curr);
        }
        // int EvalGradNext()
        // {
        //     return fatropnlp_->EvalGrad(
        //         fatropdata_->obj_scale,
        //         fatropdata_->x_next,
        //         fatropdata_->grad_next);
        // }
        double EvalObjCurr()
        {
            double res = 0.0;
            fatropnlp_->EvalObj(
                fatropdata_->obj_scale,
                fatropdata_->x_curr,
                res);
            return res;
        }
        int EvalDuInf()
        {
            return fatropnlp_->EvalDuInf(
                fatropdata_->obj_scale,
                fatropdata_->lam_curr,
                fatropdata_->grad_curr,
                fatropdata_->du_inf_curr);
        }
        inline int Initialization()
        {
            return fatropnlp_->Initializiaton(
                fatropdata_->grad_curr,
                fatropdata_->lam_calc,
                fatropdata_->delta_x);
        }
        int ComputeSD(double inertia_correction_w, double inertia_correction_c)
        {
            return fatropnlp_->ComputeSD(
                inertia_correction_w,
                inertia_correction_c,
                fatropdata_->delta_x,
                fatropdata_->lam_calc);
        }
        RefCountPtr<FatropNLP> fatropnlp_;
        RefCountPtr<FatropData> fatropdata_;
        RefCountPtr<FatropParams> fatropparams_;
        RefCountPtr<Filter> filter_;
        RefCountPtr<LineSearch> linesearch_;
        RefCountPtr<Journaller> journaller_;

    private:
        double lammax;
        double tol;
        int maxiter;
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
    };
} // namespace fatrop
#endif // FATROPALGINCLUDED