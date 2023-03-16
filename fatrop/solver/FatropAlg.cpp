#include "solver/FatropAlg.hpp"
using namespace fatrop;

FatropAlg::FatropAlg(
    const shared_ptr<FatropNLP> &fatropnlp,
    const shared_ptr<FatropData> &fatropdata,
    const shared_ptr<FatropOptions> &fatropparams,
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
    fatropnlp_->GetInitialGuess(fatropdata_->x_initial);
    fatropnlp->GetBounds(fatropdata->s_lower_orig, fatropdata->s_upper_orig);
    fatropdata->RelaxBounds();
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
    max_watchdog_steps = fatropparams_->max_watchdog_steps;
    acceptable_tol = fatropparams_->acceptable_tol;
    acceptable_iter = fatropparams_->acceptable_iter;
    warm_start_init_point = fatropparams_->warm_start_init_point;
    fatropdata_->Initialize();
    linesearch_->Initialize();
    // first_try_watchdog = fatropparams_->first_try_watchdog;
    // todo avoid reallocation when maxiter doesn't change
    // filter_ = RefCountPtr<Filter>(new Filter(maxiter + 1));
}
void FatropAlg::Reset()
{
    filter_->Reset();
    fatropdata_->Reset();
    journaller_->Reset();
    fatropnlp_->Reset();
    linesearch_->Reset();
    stats = FatropStats();
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
    // bool first_try_watchdog = this->first_try_watchdog;
    int no_watch_dog_steps_taken = 0;
    int max_watchdog_steps = this->max_watchdog_steps;
    Initialize();
    int no_conse_small_sd = false;
    int filter_reseted = 0;
    int no_no_full_steps = 0;
    int no_no_full_steps_bc_filter = 0;
    int no_acceptable_steps = 0;
    double delta_w_last_backup = 0.;
    bool restore_watchdog_step = false;
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    Reset();
    const double mu_min = tol / 10;
    double mu = mu0;
    double delta_w_last = 0.0;
    LineSearchInfo lsinfo;
    EvalJac(); // todo twice evaluation
    EvalGradCurr();
    if (warm_start_init_point)
    {
        // fatropdata_->Warm
    }
    else
    {
        int initialization_res = Initialization();
        if (initialization_res == 0 && fatropdata_->LamLinfCalc() < lammax)
        {
            cout << PRIORITY1 << "accepted lam " << endl;
            fatropdata_->AcceptInitialization();
        }
        else
        {
            cout << PRIORITY1 << "rejected lam " << endl;
            fatropdata_->lam_curr.SetConstant(0.0);
        }
    }
    fatropdata_->BoundSlacks();
    EvalCVCurr();
    fatropdata_->theta_min = fatropparams_->theta_min * MAX(1.0, fatropdata_->CVL1Curr());
    double theta_max = 1e4 * fatropdata_->CVL1Curr();
    filter_->Augment(FilterData(0, std::numeric_limits<double>::infinity(), theta_max));
    int ls = 0;
    double deltaw = 0;
    double deltac = 0.0;
    bool watch_dog_step = false;
    for (int i = 0; i < maxiter; i++)
    {
        fatropdata_->obj_curr = EvalObjCurr();
        // if (fatropdata_->LamLinfCurr() > 1e12)
        // {
        //     cout << "huge Lagrange multipliers -> set to zero" << endl;
        //     fatropdata_->lam_curr.SetConstant(0.0);
        // }
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
        if (no_no_full_steps_bc_filter >= 5)
        {
            bool reset_filter = (filter_reseted <= 5);
            if (reset_filter)
            {
                cout << PRIORITY1 << "resetted filter " << endl;
                filter_reseted++;
                filter_->Reset();
                no_no_full_steps_bc_filter = 0;
                filter_->Augment(FilterData(0, std::numeric_limits<double>::infinity(), theta_max));
            }
        }
        if (no_no_full_steps >= 10)
        {
            if (((max_watchdog_steps > 0) && !watch_dog_step))
            {
                // activate watchdog procedure
                // backup x_k
                fatropdata_->BackupCurr();
                delta_w_last_backup = delta_w_last;
                watch_dog_step = true;
                // no_no_full_steps = 0;
                no_watch_dog_steps_taken = 0;
            }
        }
        journaller_->Push();
        journaller_->PrintIterations();
        double emu = fatropdata_->EMuCurr(0.0);
        if (emu < acceptable_tol)
        {
            no_acceptable_steps++;
        }
        else
        {
            no_acceptable_steps = 0;
        }

        if (emu < tol || (no_acceptable_steps >= acceptable_iter) || ((no_conse_small_sd == 2) && (mu <= mu_min)))
        // if (emu < tol)
        {
            double total_time = blasfeo_toc(&timer);
            journaller_->PrintIterations();
            if (no_conse_small_sd == 2)
            {
                cout << "WARNING fatrop returned bc of very small search direction" << endl;
            }
            if (emu > tol && no_acceptable_steps >= acceptable_iter)
            {
                cout << "WARNING fatrop returned acceptable tolerance" << endl;
            }
            cout << "found solution :) " << endl;
            stats.eval_cv_count += linesearch_->eval_cv_count;
            stats.eval_obj_count += linesearch_->eval_obj_count;
            stats.eval_cv_time += linesearch_->eval_cv_time;
            stats.eval_obj_time += linesearch_->eval_obj_time;
            stats.time_total = total_time;
            stats.iterations_count = i;
            stats.Print();
            fatropnlp_->Finalize();
            return 0;
        }
        // update mu
        // todo make a seperate class
        while (!watch_dog_step && mu > mu_min && (fatropdata_->EMuCurr(mu) <= kappa_eta * mu || (no_conse_small_sd == 2)))
        {
            mu = MAX(mu_min, MIN(kappa_mu * mu, pow(mu, theta_mu)));
            filter_reseted = 0;
            filter_->Reset();
            no_no_full_steps_bc_filter = 0;
            filter_->Augment(FilterData(0, std::numeric_limits<double>::infinity(), theta_max));
            if (no_conse_small_sd == 2)
            {
                // cout << "small search direction" << endl;
                no_conse_small_sd = 0;
                break;
            }
            no_no_full_steps = 0;
            // the following break statement prohibits 'fast' mu updates, at leat one iteration per mu update
            // break;
        }
        // Hessian is necessary for calculating search direction
        EvalHess();
        // todo make an update class for regularization
        double deltac_candidate = delta_c_stripe * pow(mu, kappa_c);
        deltaw = 0.0;
        deltac = 0.0;
        fatropdata_->ComputeBarrierQuantities(mu);
        int regularity = -1;
        int increase_counter = 0;
        if (!restore_watchdog_step)
        {
            while (regularity != 0)
            {
                regularity = ComputeSD(deltaw, deltac, mu);
                if (deltac == 0 && regularity < 0)
                {
                    cout << PRIORITY1 << "degenerate Jacobian" << endl;
                    deltac = deltac_candidate;
                }
                if (regularity > 0) // regularization is necessary
                {
                    if (increase_counter == 0)
                    {
                        deltaw = (delta_w_last == 0.0) ? delta_w0 : MAX(delta_wmin, kappa_wmin * delta_w_last);
                    }
                    else
                    {
                        deltaw = (delta_w_last == 0.0) ? kappa_wplusem * deltaw : kappa_wplus * deltaw;
                    }
                    increase_counter++;
                }
            }
            if (deltaw > 0.)
                delta_w_last = deltaw;
        }
        fatropdata_->ComputedZ();
        // cout << "norm dzL " << Linf(fatropdata_->delta_zL) << endl;
        // cout << "norm dzU " << Linf(fatropdata_->delta_zU) << endl;
        // cout << "norm delta_s " << Linf(fatropdata_->delta_s) << endl;
        // cout << "norm delta_x " << Linf(fatropdata_->delta_x) << endl;
        // cout << "norm delta_lam " << Linf(fatropdata_->lam_calc) << endl;
        double stepsize = std::max(LinfScaled(fatropdata_->delta_x, fatropdata_->x_curr), LinfScaled(fatropdata_->delta_s, fatropdata_->s_curr));
        bool small_search_direction_curr = stepsize < 1e-14;
        lsinfo = linesearch_->FindAcceptableTrialPoint(mu, small_search_direction_curr || watch_dog_step, watch_dog_step);
        fatropdata_->RelaxBoundsVar(mu);
        fatropdata_->AdaptDualBounds(mu);
        ls = lsinfo.ls;
        if (watch_dog_step && no_watch_dog_steps_taken == 0)
        {
            fatropdata_->BackupDelta();
        }
        if (watch_dog_step)
        {
            if (ls == 1)
            {
                // accept watchdog step -- continue
                cout << PRIORITY1 << "accepted watchdog step" << endl;
                watch_dog_step = false;
            }
            else
            {
                no_watch_dog_steps_taken++;
                if (no_watch_dog_steps_taken >= max_watchdog_steps)
                {
                    // reject watchdog step -- go back to x_k
                    cout << PRIORITY1 << "rejected watchdog step" << endl;
                    it_curr.type = 'x';
                    fatropdata_->RestoreBackup();
                    // delta_w_last = delta_w_last_backup;
                    watch_dog_step = false;
                    // todo make use of delta_x_backup and delta_s_backup
                    restore_watchdog_step = true;
                    i--;
                    continue;
                };
                it_curr.type = 'w';
            }
        }

        if (ls == 0)
        {
            cout << "error: restoration phase not implemented yet" << endl;
            return 1;
        }
        if (watch_dog_step || ls == 1)
        {
            no_no_full_steps = 0;
            no_no_full_steps_bc_filter = 0;
        }
        else
        {
            ++no_no_full_steps;
            if (lsinfo.last_rejected_by_filter)
            {
                ++no_no_full_steps_bc_filter;
            }
            else
            {
                no_no_full_steps_bc_filter = 0;
            }
        }
        if (small_search_direction_curr)
        {
            it_curr.type = 's';
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
        restore_watchdog_step = false;
        // if linesearch unsuccessful -> resto phase
    }
    journaller_->PrintIterations();
    return 0;
}
int FatropAlg::EvalHess()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    int res =
        fatropnlp_->EvalHess(
            fatropdata_->obj_scale,
            fatropdata_->x_curr,
            fatropdata_->lam_curr);
    stats.eval_hess_time += blasfeo_toc(&timer);
    stats.eval_hess_count++;
    return res;
}
int FatropAlg::EvalJac()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    int res = fatropnlp_->EvalJac(
        fatropdata_->x_curr,
        fatropdata_->s_curr);
    stats.eval_jac_time += blasfeo_toc(&timer);
    stats.eval_jac_count++;
    return res;
}
inline int FatropAlg::EvalCVCurr()
{

    blasfeo_timer timer;
    blasfeo_tic(&timer);
    int res = fatropnlp_->EvalConstraintViolation(
        fatropdata_->x_curr,
        fatropdata_->s_curr,
        fatropdata_->g_curr);
    stats.eval_cv_time += blasfeo_toc(&timer);
    stats.eval_cv_count++;
    return res;
}
inline int FatropAlg::EvalGradCurr()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    int res = fatropnlp_->EvalGrad(
        fatropdata_->obj_scale,
        fatropdata_->x_curr,
        fatropdata_->grad_curr);
    stats.eval_grad_time += blasfeo_toc(&timer);
    stats.eval_grad_count++;
    return res;
}
double FatropAlg::EvalObjCurr()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    double res = 0.0;
    fatropnlp_->EvalObj(
        fatropdata_->obj_scale,
        fatropdata_->x_curr,
        res);
    stats.eval_obj_time += blasfeo_toc(&timer);
    stats.eval_obj_count++;
    return res;
}
int FatropAlg::EvalDuInf()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    fatropnlp_->EvalDuInf(
        fatropdata_->obj_scale,
        fatropdata_->lam_curr,
        fatropdata_->grad_curr,
        fatropdata_->du_inf_curr);
    fatropdata_->EvalDuInfSlacksEqs();
    stats.duinf_time += blasfeo_toc(&timer);
    return 0;
}
inline int FatropAlg::Initialization()
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    int res = fatropnlp_->Initialization_s(
        fatropdata_->s_curr);
    res = fatropnlp_->Initialization_dual(
        fatropdata_->grad_curr,
        fatropdata_->lam_calc,
        // fatropdata_->s_curr,
        fatropdata_->zL_curr,
        fatropdata_->zU_curr);
    stats.initialization_time += blasfeo_toc(&timer);
    return res;
}
int FatropAlg::ComputeSD(double inertia_correction_w, double inertia_correction_c, double mu)
{
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    int res = fatropnlp_->ComputeSD(
        inertia_correction_w,
        inertia_correction_c,
        fatropdata_->delta_x,
        fatropdata_->lam_calc,
        fatropdata_->delta_s,
        fatropdata_->sigma_total,
        fatropdata_->gradb_total);
    double el = blasfeo_toc(&timer);
    stats.compute_sd_time += el;
    return res;
}
