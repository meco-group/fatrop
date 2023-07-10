/*
 * Fatrop - A fast trajectory optimization solver
 * Copyright (C) 2022, 2023 Lander Vanroye, KU Leuven. All rights reserved.
 *
 * This file is part of Fatrop.
 *
 * Fatrop is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fatrop is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Fatrop.  If not, see <http://www.gnu.org/licenses/>. */
#include "solver/FatropData.hpp"
using namespace fatrop;
using namespace std;
template <fatrop_int size>
#define SUMMATION_ALG kahan_sum
double kahan_sum(const double *numbers)
{
    if (size == 0)
    {
        return 0.0;
    }
    double sum = numbers[0];
    double c = 0.0;
    for (fatrop_int ii = 1; ii < size; ii++)
    {
        double y = numbers[ii] - c;
        volatile double t = sum + y;
        volatile double z = t - sum;
        c = z - y;
        sum = t;
    }
    return sum;
}
template <fatrop_int size>
double kahan_sum_ld(const double *numbers)
{
    if (size == 0)
    {
        return 0.0;
    }
    long double sum = numbers[0];
    long double c = 0.0;
    for (fatrop_int ii = 1; ii < size; ii++)
    {
        long double y = numbers[ii] - c;
        volatile long double t = sum + y;
        volatile long double z = t - sum;
        c = z - y;
        sum = t;
    }
    return sum;
}
template <fatrop_int size>
double normal_sum(const double *numbers)
{
    double sum = 0.0;
    for (fatrop_int ii = 0; ii < size; ii++)
    {
        sum += numbers[ii];
    }
    return sum;
}
FatropData::FatropData(const NLPDims &nlpdims, const shared_ptr<FatropOptions> &params, const shared_ptr<FatropPrinter> &printer) : nlpdims(nlpdims),
                                                                                                                                    n_eqs(nlpdims.neqs),
                                                                                                                                    n_ineqs(nlpdims.nineqs),
                                                                                                                                    memvars(nlpdims.nvars, 12),
                                                                                                                                    memeqs(nlpdims.neqs, 12),
                                                                                                                                    memineqs(nlpdims.nineqs, 28),
                                                                                                                                    x_curr(memvars[0]),
                                                                                                                                    x_next(memvars[1]),
                                                                                                                                    x_backup(memvars[2]),
                                                                                                                                    x_initial(memvars[3]),
                                                                                                                                    delta_x(memvars[4]),
                                                                                                                                    delta_x_backup(memvars[5]),
                                                                                                                                    delta_x_backup_ls(memvars[6]),
                                                                                                                                    x_scales(memvars[7]),
                                                                                                                                    lam_curr(memeqs[0]),
                                                                                                                                    lam_next(memeqs[1]),
                                                                                                                                    lam_backup(memeqs[2]),
                                                                                                                                    lam_calc(memeqs[3]),
                                                                                                                                    lam_calc_backup(memeqs[4]),
                                                                                                                                    lam_calc_backup_ls(memeqs[5]),
                                                                                                                                    lam_scales(memeqs[6]),
                                                                                                                                    lam_init(memeqs[7]),
                                                                                                                                    g_curr(memeqs[8]),
                                                                                                                                    g_next(memeqs[9]),
                                                                                                                                    g_backup(memeqs[10]),
                                                                                                                                    g_soc(memeqs[11]),
                                                                                                                                    grad_curr(memvars[8]),
                                                                                                                                    grad_next(memvars[9]),
                                                                                                                                    grad_backup(memvars[10]),
                                                                                                                                    du_inf_curr(memvars[11]),
                                                                                                                                    du_inf_curr_s(memineqs[0]),
                                                                                                                                    s_curr(memineqs[1]),
                                                                                                                                    s_next(memineqs[2]),
                                                                                                                                    s_backup(memineqs[3]),
                                                                                                                                    delta_s(memineqs[4]),
                                                                                                                                    delta_s_backup(memineqs[5]),
                                                                                                                                    delta_s_backup_ls(memineqs[6]),
                                                                                                                                    zL_curr(memineqs[7]),
                                                                                                                                    zL_next(memineqs[8]),
                                                                                                                                    zL_backup(memineqs[9]),
                                                                                                                                    zL_init(memineqs[10]),
                                                                                                                                    zU_curr(memineqs[11]),
                                                                                                                                    zU_next(memineqs[12]),
                                                                                                                                    zU_backup(memineqs[13]),
                                                                                                                                    zU_init(memineqs[14]),
                                                                                                                                    delta_zL(memineqs[15]),
                                                                                                                                    delta_zU(memineqs[16]),
                                                                                                                                    s_lower_orig(memineqs[17]),
                                                                                                                                    s_upper_orig(memineqs[18]),
                                                                                                                                    s_lower(memineqs[19]),
                                                                                                                                    s_upper(memineqs[20]),
                                                                                                                                    sigma_L(memineqs[21]),
                                                                                                                                    sigma_U(memineqs[22]),
                                                                                                                                    sigma_total(memineqs[23]),
                                                                                                                                    gradb_L(memineqs[24]),
                                                                                                                                    gradb_U(memineqs[25]),
                                                                                                                                    gradb_plus(memineqs[26]),
                                                                                                                                    gradb_total(memineqs[27]),
                                                                                                                                    params(params),
                                                                                                                                    printer_(printer)
{
    initialize();
    params->register_option(NumericOption::lower_bounded("warm_start_mult_bound_push", "warm_start_mult_bound_push", &warm_start_mult_bound_push, 1e-2, 0.0));
    params->register_option(NumericOption::lower_bounded("smax", "smax", &smax, 100.0, 0.0));
    params->register_option(NumericOption::lower_bounded("bound_push", "kappa1", &kappa1, 1e-2, 0.0));
    params->register_option(NumericOption::lower_bounded("bound_frac", "kappa2", &kappa2, 1e-2, 0.0));
    params->register_option(NumericOption::lower_bounded("kappa_sigma", "kappa_sigma", &kappa_sigma, 1e10, 0.0));
    params->register_option(NumericOption::lower_bounded("bound_relax_factor", "bound_relax_factor", &bound_relax_factor, 1e-8, 0.0));
    params->register_option(NumericOption::lower_bounded("constr_viol_tol", "constr_viol_tol", &constr_viol_tol, 1e-4, 0.0));
}
void FatropData::initialize()
{
    kappa_d = params->kappa_d;
    n_ineqs_r = number_of_bounds();
    relax_bounds();
}
fatrop_int FatropData::reset()
{
    VEC *lower_bound_p = (VEC *)s_lower;
    VEC *upper_bound_p = (VEC *)s_upper;
    VEC *zL_p = (VEC *)zL_curr;
    VEC *zU_p = (VEC *)zU_curr;
    for (fatrop_int i = 0; i < n_ineqs; i++)
    {
        double loweri = VECEL(lower_bound_p, i);
        double upperi = VECEL(upper_bound_p, i);
        VECEL(zL_p, i) = !isinf(loweri) ? 1.0 : 0.0;
        VECEL(zU_p, i) = !isinf(upperi) ? 1.0 : 0.0;
    }
    VECSE(lam_curr.nels(), 0.0, (VEC *)lam_curr, 0);
    VECSE(s_curr.nels(), 0.0, (VEC *)s_curr, 0);
    x_curr.copy(x_initial);
    reset_caches();
    return 0;
}
fatrop_int FatropData::reset_caches()
{
    cache_curr = EvalCache();
    cache_next = EvalCache();
    return 0;
}
double FatropData::e_mu_curr(double mu)
{
    double res = 0.0;
    double z_L1 = +z_sum_curr();
    double lammean = (dual_sum_curr() + z_L1) / (n_eqs + n_ineqs_r);
    double z_mean = z_L1 / n_ineqs_r;
    double cv = constr_viol_max_curr();
    double du = dual_inf_max_curr();
    double compl_slack = eval_compl_slack(mu);
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
fatrop_int FatropData::eval_dual_inf_slack_eqs()
{
    VEC *lower_bound_p = (VEC *)s_lower;
    VEC *upper_bound_p = (VEC *)s_upper;
    VEC *lam_curr_p = (VEC *)lam_curr;
    VEC *du_inf_curr_s_p = (VEC *)du_inf_curr_s;
    VECCPSC(n_ineqs, -1.0, lam_curr_p, n_eqs - n_ineqs, du_inf_curr_s_p, 0);
    // lam_curr.print();
    VEC *zL_p = (VEC *)zL_curr;
    VEC *zU_p = (VEC *)zU_curr;
    for (fatrop_int i = 0; i < n_ineqs; i++)
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
double FatropData::eval_compl_slack(double mu)
{
    VEC *lower_bound_p = (VEC *)s_lower;
    VEC *upper_bound_p = (VEC *)s_upper;
    VEC *s_curr_p = (VEC *)s_curr;
    VEC *zL_p = (VEC *)zL_curr;
    VEC *zU_p = (VEC *)zU_curr;
    double res = 0.0;
    for (fatrop_int i = 0; i < s_curr.nels(); i++)
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
double FatropData::eval_barrier_func(double mu, VEC *s_p)
{
    VEC *lower_bound_p = (VEC *)s_lower;
    VEC *upper_bound_p = (VEC *)s_upper;
    // VEC *s_p = (VEC *)s_curr;
    double res = 0.0;
    for (fatrop_int i = 0; i < s_curr.nels(); i++)
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
double FatropData::eval_barrier_func_curr(double mu)
{
    return eval_barrier_func(mu, (VEC *)s_curr);
}
double FatropData::eval_barrier_func_backup(double mu)
{
    return eval_barrier_func(mu, (VEC *)s_backup);
}
double FatropData::eval_barrier_func_trial(double mu)
{
    return eval_barrier_func(mu, (VEC *)s_next);
}
double FatropData::eval_barrier_fo_decr(double mu, VEC *s_p, VEC *delta_s_p)
{
    VEC *lower_bound_p = (VEC *)s_lower;
    VEC *upper_bound_p = (VEC *)s_upper;
    // VEC *s_p = (VEC *)s_curr;
    // VEC *delta_s_p = (VEC *)delta_s;
    double res = 0.0;
    for (fatrop_int i = 0; i < s_curr.nels(); i++)
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
double FatropData::eval_barrier_fo_decr_curr(double mu)
{
    VEC *s_p = (VEC *)s_curr;
    VEC *delta_s_p = (VEC *)delta_s;
    return eval_barrier_fo_decr(mu, s_p, delta_s_p);
}
double FatropData::eval_barrier_fo_decr_backup(double mu)
{
    VEC *s_p = (VEC *)s_backup;
    VEC *delta_s_p = (VEC *)delta_s_backup;
    return eval_barrier_fo_decr(mu, s_p, delta_s_p);
}
fatrop_int FatropData::bound_slacks()
{
    VEC *s_lower_p = (VEC *)s_lower;
    VEC *s_upper_p = (VEC *)s_upper;
    VEC *s_curr_p = (VEC *)s_curr;
    for (fatrop_int i = 0; i < n_ineqs; i++)
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
fatrop_int FatropData::bound_z()
{
    VEC *zL_curr_p = (VEC *)zL_curr;
    VEC *zU_curr_p = (VEC *)zU_curr;
    for (fatrop_int i = 0; i < n_ineqs; i++)
    {
        VECEL(zL_curr_p, i) = MAX(VECEL(zL_curr_p, i), warm_start_mult_bound_push);
        VECEL(zU_curr_p, i) = MAX(VECEL(zU_curr_p, i), warm_start_mult_bound_push);
    }
    return 0;
}
fatrop_int FatropData::warmstart_dual()
{
    zL_curr.copy(zL_init);
    zU_curr.copy(zU_init);
    lam_curr.copy(lam_init);
    return 0;
}
fatrop_int FatropData::modify_dual_bounds(double mu)
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
    for (fatrop_int i = 0; i < n_ineqs; i++)
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
fatrop_int FatropData::accept_dual_initializiaton()
{
    lam_calc.SwapWith(lam_curr);
    cache_curr = EvalCache();
    return 0;
}
fatrop_int FatropData::update_trial_step(double alpha_primal, double alpha_dual)
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
fatrop_int FatropData::accept_trial_step()
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
fatrop_int FatropData::backup_curr()
{
    x_backup.copy(x_curr);
    s_backup.copy(s_curr);
    lam_backup.copy(lam_curr);
    zL_backup.copy(zL_curr);
    zU_backup.copy(zU_curr);
    grad_backup.copy(grad_curr);
    g_backup.copy(g_curr);
    obj_backup = obj_curr;
    return 0;
}
fatrop_int FatropData::backup_delta()
{
    delta_x_backup.copy(delta_x);
    delta_s_backup.copy(delta_s);
    lam_calc_backup.copy(lam_calc);
    return 0;
}
fatrop_int FatropData::restore_backup()
{
    x_curr.SwapWith(x_backup);
    s_curr.SwapWith(s_backup);
    lam_curr.SwapWith(lam_backup);
    zL_curr.SwapWith(zL_backup);
    zU_curr.SwapWith(zU_backup);
    grad_backup.SwapWith(grad_curr);
    g_backup.SwapWith(g_curr);
    delta_x_backup.SwapWith(delta_x);
    delta_s_backup.SwapWith(delta_s);
    lam_calc.SwapWith(lam_calc_backup);
    cache_curr = EvalCache();
    return 0;
}
double FatropData::constr_viol_max_curr()
{
    return CACHEMACRO(cache_curr.cv_linf, Linf(g_curr));
}
double FatropData::constr_viol_max_next()
{
    return CACHEMACRO(cache_next.cv_linf, Linf(g_next));
}
double FatropData::constr_viol_sum_curr()
{
    return CACHEMACRO(cache_curr.cv_l1, L1(g_curr));
}
double FatropData::constr_viol_sum_backup()
{
    return L1(g_backup);
}
double FatropData::constr_viol_sum_next()
{
    return CACHEMACRO(cache_next.cv_l1, L1(g_next));
}
double FatropData::dual_sum_curr()
{
    return CACHEMACRO(cache_curr.lam_l1, Linf(lam_curr));
}
double FatropData::dual_max_curr()
{
    return Linf(lam_curr);
}
double FatropData::dual_mean_curr()
{
    return dual_sum_curr() / nlpdims.nvars;
}
double FatropData::z_sum_curr()
{
    VEC *lower_bound_p = (VEC *)s_lower;
    VEC *upper_bound_p = (VEC *)s_upper;
    VEC *zL_p = (VEC *)zL_curr;
    VEC *zU_p = (VEC *)zU_curr;
    double res = 0.0;
    for (fatrop_int i = 0; i < n_ineqs; i++)
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
fatrop_int FatropData::number_of_bounds()
{
    VEC *lower_bound_p = (VEC *)s_lower;
    VEC *upper_bound_p = (VEC *)s_upper;
    fatrop_int res = 0;
    for (fatrop_int i = 0; i < n_ineqs; i++)
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
double FatropData::delta_dual_max()
{
    return Linf(lam_calc);
}
double FatropData::dual_inf_max_curr()
{
    return CACHEMACRO(cache_curr.du_inf_linf, MAX(Linf(du_inf_curr), Linf(du_inf_curr_s)));
}
double FatropData::fo_decr_obj_curr()
{
    return dot(grad_curr, delta_x);
}
double FatropData::fo_decr_obj_backup()
{
    return dot(grad_backup, delta_x_backup);
}
void FatropData::maximum_step_size(double &alpha_max_pr, double &alpha_max_du, double tau)
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
    for (fatrop_int i = 0; i < n_ineqs; i++)
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
void FatropData::set_bounds(const vector<double> &lowerin, const vector<double> &upperin)
{
    s_lower_orig = lowerin;
    s_upper_orig = upperin;
    relax_bounds();
}
void FatropData::relax_bounds()
{
    VECCP(n_ineqs, (VEC *)s_lower_orig, 0, (VEC *)s_lower, 0);
    VECCP(n_ineqs, (VEC *)s_upper_orig, 0, (VEC *)s_upper, 0);
    VEC *s_lower_p = (VEC *)s_lower;
    VEC *s_upper_p = (VEC *)s_upper;
    for (fatrop_int i = 0; i < n_ineqs; i++)
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
void FatropData::relax_bounds_var(double mu)
{
    double emach = 1e-16;
    VEC *s_lower_p = (VEC *)s_lower;
    VEC *s_upper_p = (VEC *)s_upper;
    VEC *s_curr_p = (VEC *)s_curr;
    for (fatrop_int i = 0; i < n_ineqs; i++)
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
                printer_->level(1) << "slacks too small " << endl;
                VECEL(s_lower_p, i) -= 1e-12 * max(1.0, std::abs(loweri));
            }
        }
        if (upper_bounded)
        {
            double s_curr_v = VECEL(s_curr_p, i);
            double dist_upper = upperi - s_curr_v;
            if (dist_upper < mu * emach)
            {
                printer_->level(1) << "slacks too small " << endl;
                VECEL(s_upper_p, i) += 1e-12 * max(1.0, std::abs(upperi));
            }
        }
    }
}
void FatropData::evaluate_barrier_quantities(double mu)
{
    VEC *s_lower_p = (VEC *)s_lower;
    VEC *s_upper_p = (VEC *)s_upper;
    VEC *s_curr_p = (VEC *)s_curr;
    VEC *zL_curr_p = (VEC *)zL_curr;
    VEC *zU_curr_p = (VEC *)zU_curr;
    VEC *sigma_L_p = (VEC *)sigma_L;
    VEC *sigma_U_p = (VEC *)sigma_U;
    VEC *gradb_L_p = (VEC *)gradb_L;
    VEC *gradb_U_p = (VEC *)gradb_U;
    VEC *gradb_plus_p = (VEC *)gradb_plus;
    VEC *lam_curr_p = (VEC *)lam_curr;
    fatrop_int eqs_offset = n_eqs - n_ineqs;
    // compute simga_LU gradb_LU and gradbplus

    for (fatrop_int i = 0; i < n_ineqs; i++)
    {
        double loweri = VECEL(s_lower_p, i);
        double upperi = VECEL(s_upper_p, i);
        bool lower_bounded = !isinf(loweri);
        bool upper_bounded = !isinf(upperi);
        if (lower_bounded)
        {
            double s_curr_v = VECEL(s_curr_p, i);
            double dist_lower = s_curr_v - loweri;
            VECEL(sigma_L_p, i) = VECEL(zL_curr_p, i) / dist_lower;
            VECEL(gradb_L_p, i) = -mu / dist_lower;
        }
        else
        {
            VECEL(sigma_L_p, i) = 0.0;
            VECEL(gradb_L_p, i) = 0.0;
        }
        if (upper_bounded)
        {
            double s_curr_v = VECEL(s_curr_p, i);
            double dist_upper = upperi - s_curr_v;
            VECEL(sigma_U_p, i) = VECEL(zU_curr_p, i) / dist_upper;
            VECEL(gradb_U_p, i) = mu / dist_upper;
        }
        else
        {
            VECEL(sigma_U_p, i) = 0.0;
            VECEL(gradb_U_p, i) = 0.0;
        }
        double grad_barrier_plusi = -VECEL(lam_curr_p, eqs_offset + i);
        if (!(lower_bounded && upper_bounded))
        {
            grad_barrier_plusi += lower_bounded ? kappa_d * mu : -kappa_d * mu;
        }
        VECEL(gradb_plus_p, i) = grad_barrier_plusi;
    }
    // total quantities
    // VECCP(n_ineqs, (VEC *)gradb_L, 0, (VEC *)gradb_total, 0);
    // AXPY(n_ineqs, 1.0, (VEC *)gradb_U, 0, (VEC *)gradb_total, 0, (VEC *)gradb_total, 0);
    // AXPY(n_ineqs, 1.0, (VEC *)gradb_plus, 0, (VEC *)gradb_total, 0, (VEC *)gradb_total, 0);
    double *gradb_total_dp = ((VEC *)gradb_total)->pa;
    double *gradb_L_dp = ((VEC *)gradb_L)->pa;
    double *gradb_U_dp = ((VEC *)gradb_U)->pa;
    double *gradb_plus_dp = ((VEC *)gradb_plus)->pa;
    for (fatrop_int i = 0; i < n_ineqs; i++)
    {
        double terms[] = {gradb_L_dp[i], gradb_U_dp[i], gradb_plus_dp[i]};
        gradb_total_dp[i] = SUMMATION_ALG<3>(terms);
    }
    VECCP(n_ineqs, (VEC *)sigma_L, 0, (VEC *)sigma_total, 0);
    AXPY(n_ineqs, 1.0, (VEC *)sigma_U, 0, (VEC *)sigma_total, 0, (VEC *)sigma_total, 0);
}
void FatropData::compute_delta_z()
{
    // delta zL
    VECMUL(n_ineqs, (VEC *)sigma_L, 0, (VEC *)delta_s, 0, (VEC *)delta_zL, 0);
    AXPBY(n_ineqs, -1.0, (VEC *)zL_curr, 0, -1.0, (VEC *)delta_zL, 0, (VEC *)delta_zL, 0);
    AXPY(n_ineqs, -1.0, (VEC *)gradb_L, 0, (VEC *)delta_zL, 0, (VEC *)delta_zL, 0);
    // delta zU
    VECMUL(n_ineqs, (VEC *)sigma_U, 0, (VEC *)delta_s, 0, (VEC *)delta_zU, 0);
    AXPBY(n_ineqs, -1.0, (VEC *)zU_curr, 0, 1.0, (VEC *)delta_zU, 0, (VEC *)delta_zU, 0);
    AXPY(n_ineqs, 1.0, (VEC *)gradb_U, 0, (VEC *)delta_zU, 0, (VEC *)delta_zU, 0);
}

void FatropData::compute_primal_dual_residu()
{
}
// void FatropData::B