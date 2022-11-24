#include "solver/FatropData.hpp"
using namespace fatrop;
FatropData::FatropData(const NLPDims &nlpdims, const shared_ptr<FatropParams> &params) : nlpdims(nlpdims),
                                                                                         n_eqs(nlpdims.neqs),
                                                                                         n_ineqs(nlpdims.nineqs),
                                                                                         memvars(nlpdims.nvars, 11),
                                                                                         memeqs(nlpdims.neqs, 8),
                                                                                         memineqs(nlpdims.nineqs, 22),
                                                                                         x_curr(memvars[0]),
                                                                                         x_next(memvars[1]),
                                                                                         x_backup(memvars[2]),
                                                                                         x_initial(memvars[3]),
                                                                                         delta_x(memvars[4]),
                                                                                         delta_x_backup(memvars[5]),
                                                                                         x_scales(memvars[6]),
                                                                                         lam_curr(memeqs[0]),
                                                                                         lam_next(memeqs[1]),
                                                                                         lam_backup(memeqs[2]),
                                                                                         lam_calc(memeqs[3]),
                                                                                         lam_scales(memeqs[4]),
                                                                                         g_curr(memeqs[5]),
                                                                                         g_next(memeqs[6]),
                                                                                         g_backup(memeqs[7]),
                                                                                         grad_curr(memvars[7]),
                                                                                         grad_next(memvars[8]),
                                                                                         grad_backup(memvars[9]),
                                                                                         du_inf_curr(memvars[10]),
                                                                                         du_inf_curr_s(memineqs[0]),
                                                                                         s_curr(memineqs[1]),
                                                                                         s_next(memineqs[2]),
                                                                                         s_backup(memineqs[3]),
                                                                                         delta_s(memineqs[4]),
                                                                                         delta_s_backup(memineqs[5]),
                                                                                         zL_curr(memineqs[6]),
                                                                                         zL_next(memineqs[7]),
                                                                                         zL_backup(memineqs[8]),
                                                                                         zU_curr(memineqs[9]),
                                                                                         zU_next(memineqs[10]),
                                                                                         zU_backup(memineqs[11]),
                                                                                         delta_zL(memineqs[12]),
                                                                                         delta_zU(memineqs[13]),
                                                                                         s_lower_orig(memineqs[14]),
                                                                                         s_upper_orig(memineqs[15]),
                                                                                         s_lower(memineqs[16]),
                                                                                         s_upper(memineqs[17]),
                                                                                         simga_L(memineqs[18]),
                                                                                         sigma_U(memineqs[19]),
                                                                                         gradb_L(memineqs[20]),
                                                                                         gradb_U(memineqs[21]),
                                                                                         gradb_plus(memineqs[22]),
                                                                                         params(params)
{
    Initialize();
}
void FatropData::Initialize()
{
    smax = params->smax;
    kappa1 = params->kappa1;
    kappa2 = params->kappa2;
    kappa_d = params->kappa_d;
    kappa_sigma = params->kappa_sigma;
    bound_relax_factor = params->bound_relax_factor;
    constr_viol_tol = params->constr_viol_tol;
    n_ineqs_r = nIneqsR();
    RelaxBounds();
}
int FatropData::Reset()
{
    VECSE(zL_curr.nels(), 1.0, (VEC *)zL_curr, 0);
    VECSE(zU_curr.nels(), 1.0, (VEC *)zU_curr, 0);
    VECSE(lam_curr.nels(), 0.0, (VEC *)lam_curr, 0);
    VECSE(s_curr.nels(), 0.0, (VEC *)s_curr, 0);
    x_curr.copy(x_initial);
    ResetCaches();
    return 0;
}
int FatropData::ResetCaches()
{
    cache_curr = EvalCache();
    cache_next = EvalCache();
    return 0;
}
double FatropData::EMuCurr(double mu)
{
    double res = 0.0;
    double z_L1 = +ZL1Curr();
    double lammean = (LamL1Curr() + z_L1) / (n_eqs + n_ineqs_r);
    double z_mean = z_L1 / n_ineqs_r;
    double cv = CVLinfCurr();
    double du = DuInfLinfCurr();
    double compl_slack = EvalCompSlackInf(mu);
    double sd = 0.0;
    double sc = 0.0;
    if (lammean > smax)
    {
        sd = lammean / smax;
        du /= sd;
    }
    if (z_mean > smax)
    {
        sc = z_mean / smax;
        compl_slack /= sc;
    }
    res = MAX(cv, MAX(du, compl_slack));
    return res;
};
int FatropData::EvalDuInfSlacksEqs()
{
    VEC *lower_bound_p = (VEC *)s_lower;
    VEC *upper_bound_p = (VEC *)s_upper;
    VEC *lam_curr_p = (VEC *)lam_curr;
    VEC *du_inf_curr_s_p = (VEC *)du_inf_curr_s;
    VECCPSC(n_ineqs, -1.0, lam_curr_p, n_eqs - n_ineqs, du_inf_curr_s_p, 0);
    // lam_curr.print();
    VEC *zL_p = (VEC *)zL_curr;
    VEC *zU_p = (VEC *)zU_curr;
    for (int i = 0; i < n_ineqs; i++)
    {
        double loweri = VECEL(lower_bound_p, i);
        double upperi = VECEL(upper_bound_p, i);
        if (!isinf(loweri))
        {
            VECEL(du_inf_curr_s_p, i) -= VECEL(zL_p, i);
        }
        if (!isinf(upperi))
        {
            VECEL(du_inf_curr_s_p, i) += VECEL(zU_p, i);
        }
    }
    return 0;
}
double FatropData::EvalCompSlackInf(double mu)
{
    VEC *lower_bound_p = (VEC *)s_lower;
    VEC *upper_bound_p = (VEC *)s_upper;
    VEC *s_curr_p = (VEC *)s_curr;
    VEC *zL_p = (VEC *)zL_curr;
    VEC *zU_p = (VEC *)zU_curr;
    double res = 0.0;
    for (int i = 0; i < s_curr.nels(); i++)
    {
        double loweri = VECEL(lower_bound_p, i);
        double upperi = VECEL(upper_bound_p, i);
        double si = VECEL(s_curr_p, i);
        if (!isinf(loweri))
        {
            double dist = si - loweri;
            res = MAX(res, dist * VECEL(zL_p, i) - mu);
        }
        if (!isinf(upperi))
        {
            double dist = upperi - si;
            res = MAX(res, dist * VECEL(zU_p, i) - mu);
        }
    }
    return res;
}
double FatropData::EvalBarrier(double mu, VEC *s_p)
{
    VEC *lower_bound_p = (VEC *)s_lower;
    VEC *upper_bound_p = (VEC *)s_upper;
    // VEC *s_p = (VEC *)s_curr;
    double res = 0.0;
    for (int i = 0; i < s_curr.nels(); i++)
    {
        double loweri = VECEL(lower_bound_p, i);
        double upperi = VECEL(upper_bound_p, i);
        double si = VECEL(s_p, i);
        bool lower_bounded = !isinf(loweri);
        bool upper_bounded = !isinf(upperi);
        bool one_sided = !(lower_bounded && upper_bounded);
        if (lower_bounded)
        {
            double dist = si - loweri;
            res += -mu * log(dist);
            if (one_sided)
                res += kappa_d * mu * dist;
        }
        if (upper_bounded)
        {
            double dist = upperi - si;
            res += -mu * log(dist);
            if (one_sided)
                res += kappa_d * mu * dist;
        }
    }
    return res;
}
double FatropData::EvalBarrierCurr(double mu)
{
    return EvalBarrier(mu, (VEC *)s_curr);
}
double FatropData::EvalBarrierBackup(double mu)
{
    return EvalBarrier(mu, (VEC *)s_backup);
}
double FatropData::EvalBarrierNext(double mu)
{
    return EvalBarrier(mu, (VEC *)s_next);
}
double FatropData::EvalBarrierLinDecr(double mu, VEC *s_p, VEC *delta_s_p)
{
    VEC *lower_bound_p = (VEC *)s_lower;
    VEC *upper_bound_p = (VEC *)s_upper;
    // VEC *s_p = (VEC *)s_curr;
    // VEC *delta_s_p = (VEC *)delta_s;
    double res = 0.0;
    for (int i = 0; i < s_curr.nels(); i++)
    {
        double loweri = VECEL(lower_bound_p, i);
        double upperi = VECEL(upper_bound_p, i);
        double si = VECEL(s_p, i);
        double delta_si = VECEL(delta_s_p, i);
        bool lower_bounded = !isinf(loweri);
        bool upper_bounded = !isinf(upperi);
        bool one_sided = !(lower_bounded && upper_bounded);
        if (lower_bounded)
        {
            double dist = si - loweri;
            res += -mu * delta_si / dist;
            if (one_sided)
                res += kappa_d * mu * delta_si;
        }
        if (upper_bounded)
        {
            double dist = upperi - si;
            res += mu * delta_si / dist;
            if (one_sided)
                res -= kappa_d * mu * delta_si;
        }
    }
    return res;
}
double FatropData::EvalBarrierLinDecrCurr(double mu)
{
    VEC *s_p = (VEC *)s_curr;
    VEC *delta_s_p = (VEC *)delta_s;
    return EvalBarrierLinDecr(mu, s_p, delta_s_p);
}
double FatropData::EvalBarrierLinDecrBackup(double mu)
{
    VEC *s_p = (VEC *)s_backup;
    VEC *delta_s_p = (VEC *)delta_s_backup;
    return EvalBarrierLinDecr(mu, s_p, delta_s_p);
}
int FatropData::BoundSlacks()
{
    VEC *s_lower_p = (VEC *)s_lower;
    VEC *s_upper_p = (VEC *)s_upper;
    VEC *s_curr_p = (VEC *)s_curr;
    for (int i = 0; i < n_ineqs; i++)
    {
        double loweri = VECEL(s_lower_p, i);
        double upperi = VECEL(s_upper_p, i);
        bool lower_bounded = !isinf(loweri);
        bool upper_bounded = !isinf(upperi);
        bool two_sided = lower_bounded && upper_bounded;
        if (two_sided)
        {
            double pL = MIN(kappa1 * MAX(1.0, abs(loweri)), kappa2 * (upperi - loweri));
            double pR = MIN(kappa1 * MAX(1.0, abs(upperi)), kappa2 * (upperi - loweri));
            DBGASSERT((pL > 0) && (pR > 0));
            VECEL(s_curr_p, i) = MIN(MAX(VECEL(s_curr_p, i), loweri + pL), upperi - pR);
            DBGASSERT((VECEL(s_curr_p, i) > loweri) || (VECEL(s_curr_p, i) < upperi));
        }
        else if (lower_bounded)
        {
            VECEL(s_curr_p, i) = MAX(VECEL(s_curr_p, i), loweri + kappa1 * MAX(1.0, abs(loweri)));
        }
        else if (upper_bounded)
        {
            VECEL(s_curr_p, i) = MIN(VECEL(s_curr_p, i), upperi - kappa1 * MAX(1.0, abs(upperi)));
        }
    }
    return 0;
}
int FatropData::AdaptDualBounds(double mu)
{
    VEC *s_lower_p = (VEC *)s_lower;
    VEC *s_upper_p = (VEC *)s_upper;
    VEC *s_curr_p = (VEC *)s_curr;
    // VEC *s_curr_p = (VEC *)s_curr;
    // VEC *delta_s_p = (VEC *)delta_s;
    // VEC * zL_curr_p = (VEC* )zL_curr;
    // VEC * zU_curr_p = (VEC* )zU_curr;
    VEC *zL_curr_p = (VEC *)zL_curr;
    VEC *zU_curr_p = (VEC *)zU_curr;
    double kappa_sigma = this->kappa_sigma;
    for (int i = 0; i < n_ineqs; i++)
    {
        double loweri = VECEL(s_lower_p, i);
        double upperi = VECEL(s_upper_p, i);
        bool lower_bounded = !isinf(loweri);
        bool upper_bounded = !isinf(upperi);
        if (lower_bounded)
        {
            double s_curr_v = VECEL(s_curr_p, i);
            double dist_lower = s_curr_v - VECEL(s_lower_p, i);
            double zL_curr_v = VECEL(zL_curr_p, i);
            VECEL(zL_curr_p, i) = MAX(MIN(zL_curr_v, kappa_sigma * mu / dist_lower), mu / (kappa_sigma * dist_lower));
        }
        if (upper_bounded)
        {
            double s_curr_v = VECEL(s_curr_p, i);
            double dist_upper = VECEL(s_upper_p, i) - s_curr_v;
            double zU_curr_v = VECEL(zU_curr_p, i);
            VECEL(zU_curr_p, i) = MAX(MIN(zU_curr_v, kappa_sigma * mu / dist_upper), mu / (kappa_sigma * dist_upper));
        }
    }
    return 0;
}
int FatropData::AcceptInitialization()
{
    lam_calc.SwapWith(lam_curr);
    cache_curr = EvalCache();
    return 0;
}
int FatropData::TryStep(double alpha_primal, double alpha_dual)
{
    axpy(alpha_primal, delta_x, x_curr, x_next);
    axpy(alpha_primal, delta_s, s_curr, s_next);
    axpy(alpha_dual, delta_zL, zL_curr, zL_next);
    axpy(alpha_dual, delta_zU, zU_curr, zU_next);
    axpy(alpha_primal, lam_calc, lam_curr, lam_next);
    // axpby(alpha_primal, lam_calc, 1.0 - alpha_primal, lam_curr, lam_next);
    // reset evaluation flags
    cache_next = EvalCache();
    return 0;
}
int FatropData::TakeStep()
{
    // TODO make a struct which containts vectors associated with curr <-> next
    x_curr.SwapWith(x_next);
    s_curr.SwapWith(s_next);
    lam_curr.SwapWith(lam_next);
    zL_curr.SwapWith(zL_next);
    zU_curr.SwapWith(zU_next);
    grad_curr.SwapWith(grad_next);
    g_curr.SwapWith(g_next);
    cache_curr = cache_next;
    return 0;
}
int FatropData::BackupCurr()
{
    x_backup.copy(x_curr);
    s_backup.copy(s_curr);
    lam_backup.copy(lam_curr);
    zL_backup.copy(zL_curr);
    zU_backup.copy(zU_curr);
    grad_backup.copy(grad_curr);
    g_backup.copy(g_curr);
    delta_x_backup.copy(delta_x);
    delta_s_backup.copy(delta_s);
    obj_backup = obj_curr;
    return 0;
}
int FatropData::RestoreBackup()
{
    x_curr.SwapWith(x_backup);
    s_curr.SwapWith(s_backup);
    lam_curr.SwapWith(lam_backup);
    zL_curr.SwapWith(zL_backup);
    zU_curr.SwapWith(zU_backup);
    grad_backup.SwapWith(grad_curr);
    g_backup.SwapWith(g_curr);
    cache_curr = EvalCache();
    return 0;
}
double FatropData::CVLinfCurr()
{
    return CACHEMACRO(cache_curr.cv_linf, Linf(g_curr));
}
double FatropData::CVLinfNext()
{
    return CACHEMACRO(cache_next.cv_linf, Linf(g_next));
}
double FatropData::CVL1Curr()
{
    return CACHEMACRO(cache_curr.cv_l1, L1(g_curr));
}
double FatropData::CVL1Backup()
{
    return L1(g_backup);
}
double FatropData::CVL1Next()
{
    return CACHEMACRO(cache_next.cv_l1, L1(g_next));
}
double FatropData::LamL1Curr()
{
    return CACHEMACRO(cache_curr.lam_l1, Linf(lam_curr));
}
double FatropData::LamLinfCurr()
{
    return Linf(lam_curr);
}
double FatropData::LamMeanCurr()
{
    return LamL1Curr() / nlpdims.nvars;
}
double FatropData::ZL1Curr()
{
    VEC *lower_bound_p = (VEC *)s_lower;
    VEC *upper_bound_p = (VEC *)s_upper;
    VEC *zL_p = (VEC *)zL_curr;
    VEC *zU_p = (VEC *)zU_curr;
    double res = 0.0;
    for (int i = 0; i < n_ineqs; i++)
    {
        double loweri = VECEL(lower_bound_p, i);
        double upperi = VECEL(upper_bound_p, i);
        if (!isinf(loweri))
        {
            res += abs(VECEL(zL_p, i));
        }
        if (!isinf(upperi))
        {
            res += abs(VECEL(zU_p, i));
        }
    }
    return res;
}
int FatropData::nIneqsR()
{
    VEC *lower_bound_p = (VEC *)s_lower;
    VEC *upper_bound_p = (VEC *)s_upper;
    int res = 0;
    for (int i = 0; i < n_ineqs; i++)
    {
        double loweri = VECEL(lower_bound_p, i);
        double upperi = VECEL(upper_bound_p, i);
        if (!isinf(loweri))
        {
            res++;
        }
        if (!isinf(upperi))
        {
            res++;
        }
    }
    return res;
}
double FatropData::LamLinfCalc()
{
    return Linf(lam_calc);
}
double FatropData::DuInfLinfCurr()
{
    return CACHEMACRO(cache_curr.du_inf_linf, MAX(Linf(du_inf_curr), Linf(du_inf_curr_s)));
}
double FatropData::LinDecrCurr()
{
    return dot(grad_curr, delta_x);
}
double FatropData::LinDecrBackup()
{
    return dot(grad_backup, delta_x_backup);
}
void FatropData::AlphaMax(double &alpha_max_pr, double &alpha_max_du, double tau)
{
    alpha_max_pr = 1.0;
    alpha_max_du = 1.0;
    VEC *s_lower_p = (VEC *)s_lower;
    VEC *s_upper_p = (VEC *)s_upper;
    VEC *delta_s_p = (VEC *)delta_s;
    VEC *s_curr_p = (VEC *)s_curr;
    VEC *zL_curr_p = (VEC *)zL_curr;
    VEC *zU_curr_p = (VEC *)zU_curr;
    VEC *delta_zL_p = (VEC *)delta_zL;
    VEC *delta_zU_p = (VEC *)delta_zU;
    for (int i = 0; i < n_ineqs; i++)
    {
        if (!isinf(VECEL(s_lower_p, i)))
        {
            double delta_s_i = VECEL(delta_s_p, i);
            double delta_Z_i = VECEL(delta_zL_p, i);
            // primal
            alpha_max_pr = delta_s_i < 0 ? MIN(alpha_max_pr, -tau * (VECEL(s_curr_p, i) - VECEL(s_lower_p, i)) / delta_s_i) : alpha_max_pr;
            // dual
            alpha_max_du = delta_Z_i < 0 ? MIN(alpha_max_du, -tau * (VECEL(zL_curr_p, i)) / delta_Z_i) : alpha_max_du;
        }
        if (!isinf(VECEL(s_upper_p, i)))
        {
            double delta_s_i = VECEL(delta_s_p, i);
            double delta_Z_i = VECEL(delta_zU_p, i);
            // primal
            alpha_max_pr = delta_s_i > 0 ? MIN(alpha_max_pr, tau * (VECEL(s_upper_p, i) - VECEL(s_curr_p, i)) / delta_s_i) : alpha_max_pr;
            // dual
            alpha_max_du = delta_Z_i < 0 ? MIN(alpha_max_du, -tau * (VECEL(zU_curr_p, i)) / delta_Z_i) : alpha_max_du;
        }
    }
    return;
}
void FatropData::SetBounds(const vector<double> &lowerin, const vector<double> &upperin)
{
    s_lower_orig = lowerin;
    s_upper_orig = upperin;
    RelaxBounds();
}
void FatropData::RelaxBounds()
{
    VECCP(n_ineqs, (VEC *)s_lower_orig, 0, (VEC *)s_lower, 0);
    VECCP(n_ineqs, (VEC *)s_upper_orig, 0, (VEC *)s_upper, 0);
    VEC *s_lower_p = (VEC *)s_lower;
    VEC *s_upper_p = (VEC *)s_upper;
    for (int i = 0; i < n_ineqs; i++)
    {
        double lower_i = VECEL(s_lower_p, i);
        if (!isinf(lower_i))
        {
            VECEL(s_lower_p, i) = lower_i - MIN(constr_viol_tol, bound_relax_factor * MAX(1.00, abs(lower_i)));
        }
        double upper_i = VECEL(s_upper_p, i);
        if (!isinf(upper_i))
        {
            VECEL(s_upper_p, i) = upper_i + MIN(constr_viol_tol, bound_relax_factor * MAX(1.00, abs(upper_i)));
        }
    }
}
void FatropData::RelaxBoundsVar(double mu)
{
    double emach = 1e-16;
    VEC *s_lower_p = (VEC *)s_lower;
    VEC *s_upper_p = (VEC *)s_upper;
    VEC *s_curr_p = (VEC *)s_curr;
    for (int i = 0; i < n_ineqs; i++)
    {
        double loweri = VECEL(s_lower_p, i);
        double upperi = VECEL(s_upper_p, i);
        bool lower_bounded = !isinf(loweri);
        bool upper_bounded = !isinf(upperi);
        if (lower_bounded)
        {
            double s_curr_v = VECEL(s_curr_p, i);
            double dist_lower = s_curr_v - loweri;
            if (dist_lower < mu * emach)
            {
                cout << "slacks too small " << endl;
                VECEL(s_lower_p, i) -= 1e-12 * std::max(1.0, std::abs(loweri));
            }
        }
        if (upper_bounded)
        {
            double s_curr_v = VECEL(s_curr_p, i);
            double dist_upper = upperi - s_curr_v;
            if (dist_upper < mu * emach)
            {
                cout << "slacks too small " << endl;
                VECEL(s_upper_p, i) += 1e-12 * std::max(1.0, std::abs(upperi));
            }
        }
    }
}
// void FatropData::B