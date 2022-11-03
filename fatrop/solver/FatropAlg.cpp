#include "solver/FatropAlg.hpp"
using namespace fatrop;

FatropAlg::FatropAlg(
    const shared_ptr<FatropNLP> &fatropnlp,
    const shared_ptr<FatropData> &fatropdata,
    const shared_ptr<FatropParams> &fatropparams,
    const shared_ptr<Filter> &filter,
    const shared_ptr<LineSearch> &linesearch,
    const shared_ptr<Journaller> &journaller)
    : fatropnlp_(fatropnlp),
      fatropdata_(fatropdata),
      fatropparams_(fatropparams),
      filter_(filter),
      linesearch_(linesearch),
      journaller_(journaller)
{
    Initialize();
}
void FatropAlg::Initialize()
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
    kappa_d = fatropparams_->kappa_d;
    // todo avoid reallocation when maxiter doesn't change
    // filter_ = RefCountPtr<Filter>(new Filter(maxiter + 1));
}
void FatropAlg::Reset()
{
    filter_->Reset();
    fatropdata_->Reset();
    journaller_->Reset();
    fatropnlp_->Reset();
    sd_time = 0.0;
    init_time = 0.0;
    total_time = 0.0;
    hess_time = 0.0;
}
void FatropAlg::SetBounds(const vector<double> &lower, const vector<double> &upper)
{
    fatropdata_->s_lower_orig = lower;
    fatropdata_->s_upper_orig = upper;
};
void FatropAlg::SetInitial(const vector<double> &initial)
{
    fatropdata_->x_initial = initial;
};
void FatropAlg::GetSolution(vector<double> &sol)
{
    fatropdata_->x_curr.copyto(sol);
};
int FatropAlg::Optimize()
{
    Initialize();
    int no_conse_small_sd = false;
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    Reset();
    // optimization algorithm parameters
    const double mu_min = tol / 10;
    // optimization variables
    double mu = mu0;
    double delta_w_last = 0.0;
    // initialization
#ifdef ENABLE_MULTITHREADING // TODO control cores to which threads are assigned and take into account hyperthreading in this.
    // TODO check if it is more interesting to make worker a member variable, instead of a local variable in Optimize
    // TODO parallelize functions themselves, openmp might be a good framework, to be futher investigated
    Worker worker([&]()
                  { this->EvalHess(); });
    EvalJac(); // todo twice evaluation
    EvalGradCurr();
#else
    EvalJac(); // todo twice evaluation
    EvalGradCurr();
#endif
    int initialization_res = Initialization();
    // cout << " initializaiton res " << initialization_res << endl;
    fatropdata_->BoundSlacks();
    if (initialization_res == 0 && fatropdata_->LamLinfCalc() < lammax)
    {
        cout << "accepted lam " << endl;
        fatropdata_->AcceptInitialization();
    }
    else
    {
        cout << "rejected lam " << endl;
        fatropdata_->lam_curr.SetConstant(0.0);
    }
    EvalCVCurr();
    fatropdata_->theta_min = fatropparams_->theta_min * MAX(1.0, fatropdata_->CVL1Curr());
    int ls = 0;
    double deltaw = 0;
    double deltac = 0.0;
    for (int i = 0; i < maxiter; i++)
    {
        // cout << "iteration " << i << endl;
        // cout << "x_curr " << endl;
        // fatropdata_-> x_curr.print();
        // cout << "s_curr " << endl;
        // fatropdata_-> s_curr.print(); //
        // cout << "lam_curr " << endl;
        // fatropdata_-> lam_curr.print();
        // cout << "zL_curr " << endl;
        // fatropdata_-> zL_curr.print();
        // cout << "zU_curr " << endl;
        // fatropdata_-> zU_curr.print();
        fatropdata_->obj_curr = EvalObjCurr();
#ifdef ENABLE_MULTITHREADING
        worker.eval();
        EvalJac();
        EvalGradCurr();
#else
        if (fatropdata_->LamLinfCurr() > 1e8)
        {
            cout << "huge Lagrange multipliers -> set to zero" << endl;
            fatropdata_->lam_curr.SetConstant(0.0);
        }
        fatropdata_->RelaxBoundsVar(mu);
        EvalJac();      // needed for dual inf
        EvalGradCurr(); // needed for dual inf
#endif
        EvalDuInf();
        IterationData &it_curr = journaller_->it_curr;
        it_curr.iter = i;
        it_curr.mu = mu;
        it_curr.objective = EvalObjCurr();
        it_curr.constraint_violation = fatropdata_->CVLinfCurr();
        it_curr.du_inf = fatropdata_->DuInfLinfCurr();
        it_curr.ls = ls;
        it_curr.reg = deltaw;
        // cout << Linf(fatropdata_->zL_curr)<< endl;
        // cout << Linf(fatropdata_->zU_curr)<< endl;
        journaller_->Push();
        journaller_->PrintIterations();
        double emu = fatropdata_->EMuCurr(0.0);
        if (emu < tol || ((no_conse_small_sd == 2) && (mu <= mu_min)))
        {
            total_time = blasfeo_toc(&timer);
            journaller_->PrintIterations();
            if (no_conse_small_sd == 2)
            {
                cout << "WARNING fatrop returned bc of very small search direction" << endl;
            }
            cout << "found solution :) " << endl;
            cout << "riccati time: " << sd_time << endl;
            cout << "init time: " << init_time << endl;
            cout << "hess time " << hess_time << endl;
            // cout << "fe time nlp_f: "
            //      << "to be added" << endl;
            // cout << "fe time nlp_g: "
            //      << "to be added" << endl;
            // cout << "fe time nlp_grad_f: "
            //      << "to be added" << endl;
            // cout << "fe time nlp_hess_l: "
            //      << "to be added" << endl;
            // cout << "fe time nlp_jac_g: "
            //      << "to be added" << endl;
            fatropnlp_->Finalize();
            // cout << "rest time: " << total_time - sd_time - init_time << endl;
            cout << "el time total: " << total_time << endl;
#ifdef ENABLE_MULTITHREADING
            worker.wait();
#endif
            return 0;
        }
        // update mu
        // todo make a seperate class
        while (mu > mu_min && (fatropdata_->EMuCurr(mu) <= kappa_eta * mu || (no_conse_small_sd == 2)))
        {
            mu = MAX(mu_min, MIN(kappa_mu * mu, pow(mu, theta_mu)));
            filter_->Reset();
            if (no_conse_small_sd == 2)
            {
                // cout << "small search direction" << endl;
                no_conse_small_sd = 0;
                break;
            }
        }
        // Hessian is necessary for calculating search direction
#ifdef ENABLE_MULTITHREADING
        worker.wait();
#else
        EvalHess();
#endif
        // todo make an update class for regularization
        double deltac_candidate = delta_c_stripe * pow(mu, kappa_c);
        deltaw = 0.0;
        // deltaw = 1e-8*mu;
        deltac = 0.0;
        int regularity = ComputeSD(deltaw, deltac, mu);
        if (regularity < 0)
        {
            deltac = deltac_candidate;
            regularity = ComputeSD(deltaw, deltac, mu);
            cout << "Jac degenerate" << endl;
        }
        int increase_counter = 0;
        if (regularity > 0) // regularization is necessary
        {
            deltaw = (delta_w_last == 0.0) ? delta_w0 : MAX(delta_wmin, kappa_wmin * delta_w_last);
            regularity = ComputeSD(deltaw, deltac, mu);
            if ((deltac == 0.0) && (regularity < 0))
            {
                deltac = deltac_candidate;
                regularity = ComputeSD(deltaw, deltac, mu);
                cout << "Jac degenerate" << endl;
            }
            while (regularity > 0)
            {
                increase_counter++;
                deltaw = (delta_w_last == 0.0) ? kappa_wplusem * deltaw : kappa_wplus * deltaw;
                regularity = ComputeSD(deltaw, deltac, mu);
                if ((deltac == 0.0) && (regularity < 0))
                {
                    deltac = deltac_candidate;
                    regularity = ComputeSD(deltaw, deltac, mu);
                    cout << "Jac degenerate" << endl;
                }
            }
            delta_w_last = deltaw;
        }
        double stepsize = std::max(Linf(fatropdata_->delta_x), Linf(fatropdata_->delta_s));
        // cout << "step size " <<stepsize  << endl;
        bool small_search_direction_curr = stepsize < 1e-6;
        // cout << "regularization  " << (deltaw) << endl;
        // cout << "step size " << Linf(fatropdata_->delta_x) << endl;
        ls = linesearch_->FindAcceptableTrialPoint(mu, small_search_direction_curr);
        if (ls == 0)
        {
            cout << "error: restoration phase not implemented yet" << endl;
            return 1;
        }
        fatropdata_->AdaptDualBounds(mu);
        if (small_search_direction_curr)
        {
            no_conse_small_sd++;
            if (!no_conse_small_sd)
            {
                // take the full step
                cout << "full step small sd " << endl;
            }
        }
        else
        {
            no_conse_small_sd = 0;
        }
        // if linesearch unsuccessful -> resto phase
    }
    journaller_->PrintIterations();
    return 0;
}
inline int FatropAlg::EvalHess()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    int res =
        fatropnlp_->EvalHess(
            fatropdata_->obj_scale,
            fatropdata_->x_curr,
            fatropdata_->lam_curr);
    hess_time += blasfeo_toc(&timer);
    return res;
}
inline int FatropAlg::EvalJac()
{
    int res = fatropnlp_->EvalJac(
        fatropdata_->x_curr,
        fatropdata_->s_curr);
    return res;
}
inline int FatropAlg::EvalCVCurr()
{
    int res = fatropnlp_->EvalConstraintViolation(
        fatropdata_->x_curr,
        fatropdata_->s_curr,
        fatropdata_->g_curr);
    return res;
}
inline int FatropAlg::EvalGradCurr()
{
    int res = fatropnlp_->EvalGrad(
        fatropdata_->obj_scale,
        fatropdata_->x_curr,
        fatropdata_->grad_curr);
    return res;
}
double FatropAlg::EvalObjCurr()
{
    double res = 0.0;
    fatropnlp_->EvalObj(
        fatropdata_->obj_scale,
        fatropdata_->x_curr,
        res);
    return res;
}
int FatropAlg::EvalDuInf()
{
    fatropnlp_->EvalDuInf(
        fatropdata_->obj_scale,
        fatropdata_->lam_curr,
        fatropdata_->grad_curr,
        fatropdata_->du_inf_curr);
    fatropdata_->EvalDuInfSlacksEqs();
    return 0;
}
inline int FatropAlg::Initialization()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    int res = fatropnlp_->Initialization(
        fatropdata_->grad_curr,
        fatropdata_->lam_calc,
        fatropdata_->delta_x,
        fatropdata_->delta_s,
        fatropdata_->s_curr,
        fatropdata_->zL_curr,
        fatropdata_->zU_curr,
        fatropdata_->s_lower,
        fatropdata_->s_upper);
    init_time += blasfeo_toc(&timer);
    return res;
}
int FatropAlg::ComputeSD(double inertia_correction_w, double inertia_correction_c, double mu)
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    int res = fatropnlp_->ComputeSD(
        inertia_correction_w,
        inertia_correction_c,
        mu,
        kappa_d,
        fatropdata_->delta_x,
        fatropdata_->lam_calc,
        fatropdata_->lam_curr,
        fatropdata_->s_curr,
        fatropdata_->zL_curr,
        fatropdata_->zU_curr,
        fatropdata_->delta_zL,
        fatropdata_->delta_zU,
        fatropdata_->s_lower,
        fatropdata_->s_upper,
        fatropdata_->delta_s);
    sd_time += blasfeo_toc(&timer);
    return res;
}
