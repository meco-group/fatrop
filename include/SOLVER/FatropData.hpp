// solver data
#ifndef FATROPDATAINCLUDED
#define FATROPDATAINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "TEMPLATES/NLPAlg.hpp"
#include "AUX/Common.hpp"
#include "FatropParams.hpp"
#include <cmath>
using namespace std;
namespace fatrop
{
#define CACHEMACRO(instance, val) instance.evaluated ? instance.value : instance.SetValue(val)
    struct FatropData : public RefCountedObj
    {
        FatropData(const NLPDims &nlpdims, const RefCountPtr<FatropParams> &params) : nlpdims(nlpdims),
                                                                                      n_eqs(nlpdims.neqs),
                                                                                      n_ineqs(nlpdims.nineqs),
                                                                                      memvars(nlpdims.nvars, 7),
                                                                                      memeqs(nlpdims.neqs, 6),
                                                                                      memineqs(nlpdims.nineqs, 12),
                                                                                      x_curr(memvars[0]),
                                                                                      x_next(memvars[1]),
                                                                                      delta_x(memvars[2]),
                                                                                      x_scales(memvars[3]),
                                                                                      lam_curr(memeqs[0]),
                                                                                      lam_next(memeqs[1]),
                                                                                      lam_calc(memeqs[2]),
                                                                                      lam_scales(memeqs[3]),
                                                                                      g_curr(memeqs[4]),
                                                                                      g_next(memeqs[5]),
                                                                                      grad_curr(memvars[4]),
                                                                                      grad_next(memvars[5]),
                                                                                      du_inf_curr(memvars[6]),
                                                                                      du_inf_curr_s(memineqs[0]),
                                                                                      s_curr(memineqs[1]),
                                                                                      s_next(memineqs[2]),
                                                                                      delta_s(memineqs[3]),
                                                                                      zL_curr(memineqs[4]),
                                                                                      zL_next(memineqs[5]),
                                                                                      zU_curr(memineqs[6]),
                                                                                      zU_next(memineqs[7]),
                                                                                      delta_zL(memineqs[8]),
                                                                                      delta_zU(memineqs[9]),
                                                                                      s_lower(memineqs[10]),
                                                                                      s_upper(memineqs[11]),
                                                                                      params(params)
        {
            Initialize();
        }
        void Initialize()
        {
            smax = params->smax;
            kappa1 = params->kappa1;
            kappa2 = params->kappa2;
            kappa_d = params->kappa_d;
            n_ineqs_r = nIneqsR();
        }
        int Reset()
        {
            VECSE(zL_curr.nels(), 1.0, (VEC *)zL_curr, 0);
            VECSE(zU_curr.nels(), 1.0, (VEC *)zU_curr, 0);
            VECSE(lam_curr.nels(), 0.0, (VEC *)lam_curr, 0);
            VECSE(s_curr.nels(), 0.0, (VEC *)s_curr, 0);
            cache_curr = EvalCache();
            cache_next = EvalCache();
            return 0;
        }
        double EMuCurr(double mu)
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
        int EvalDuInfSlacksEqs()
        {
            VEC *lower_bound_p = (VEC *)s_lower;
            VEC *upper_bound_p = (VEC *)s_upper;
            VEC *lam_curr_p = (VEC *)lam_curr;
            VEC *du_inf_curr_s_p = (VEC *)du_inf_curr_s;
            VECCPSC(n_ineqs, -1.0, lam_curr_p, n_eqs, du_inf_curr_s_p, 0);
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
                    VECEL(du_inf_curr_s_p, i) -= VECEL(zU_p, i);
                }
            }
            return 0;
        }
        double EvalCompSlackInf(double mu)
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
        double EvalBarrier(double mu, VEC *s_p)
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
        double EvalBarrierCurr(double mu)
        {
            return EvalBarrier(mu, (VEC *)s_curr);
        }
        double EvalBarrierNext(double mu)
        {
            return EvalBarrier(mu, (VEC *)s_next);
        }
        double EvalBarrierLinDecr(double mu)
        {
            VEC *lower_bound_p = (VEC *)s_lower;
            VEC *upper_bound_p = (VEC *)s_upper;
            VEC *s_p = (VEC *)s_curr;
            VEC *delta_s_p = (VEC *)delta_s;
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
                        res += kappa_d * mu*delta_si;
                }
                if (upper_bounded)
                {
                    double dist = upperi - si;
                    res += mu * delta_si / dist;
                    if (one_sided)
                        res -= kappa_d * mu*delta_si;
                }
            }
            return res;
        }
        int BoundSlacks()
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
                    // project
                    VECEL(s_curr_p, i) = MIN(MAX(VECEL(s_curr_p, i), loweri + pL), upperi - pR);
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
        int AcceptInitialization()
        {
            lam_calc.SwapWith(lam_curr);
            cache_curr = EvalCache();
            return 0;
        }
        int TryStep(double alpha_primal, double alpha_dual)
        {
            axpy(alpha_primal, delta_x, x_curr, x_next);
            axpy(alpha_primal, delta_s, s_curr, s_next);
            axpy(alpha_dual, delta_zL, zL_curr, zL_next);
            axpy(alpha_dual, delta_zU, zU_curr, zU_next);
            axpby(alpha_primal, lam_calc, 1.0 - alpha_primal, lam_curr, lam_next);
            // reset evaluation flags
            cache_next = EvalCache();
            return 0;
        }
        int TakeStep()
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
        double CVLinfCurr()
        {
            return CACHEMACRO(cache_curr.cv_linf, Linf(g_curr));
        }
        double CVLinfNext()
        {
            return CACHEMACRO(cache_next.cv_linf, Linf(g_next));
        }
        double CVL1Curr()
        {
            return CACHEMACRO(cache_curr.cv_l1, L1(g_curr));
        }
        double CVL1Next()
        {
            return CACHEMACRO(cache_next.cv_l1, L1(g_next));
        }
        double LamL1Curr()
        {
            return CACHEMACRO(cache_curr.lam_l1, Linf(lam_curr));
        }
        double LamLinfCurr()
        {
            return Linf(lam_curr);
        }
        double LamMeanCurr()
        {
            return LamL1Curr() / nlpdims.nvars;
        }
        double ZL1Curr()
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
        int nIneqsR()
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
        double LamLinfCalc()
        {
            return Linf(lam_calc);
        }
        double DuInfLinfCurr()
        {
            return CACHEMACRO(cache_curr.du_inf_linf, Linf(du_inf_curr));
        }
        double LinDecrCurr()
        {
            return dot(grad_curr, delta_x);
        }
        void AlphaMax(double &alpha_max_pr, double &alpha_max_du, double tau)
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

        const NLPDims nlpdims;
        double obj_scale = 1.0;
        int n_eqs;
        int n_ineqs;
        int n_ineqs_r = 0;
        FatropMemoryVecBF memvars;
        FatropMemoryVecBF memeqs;
        FatropMemoryVecBF memineqs;
        FatropVecBF x_curr;
        FatropVecBF x_next;
        FatropVecBF delta_x;
        FatropVecBF x_scales;
        FatropVecBF lam_curr;
        FatropVecBF lam_next;
        FatropVecBF lam_calc;
        FatropVecBF lam_scales;
        FatropVecBF g_curr;
        FatropVecBF g_next;
        FatropVecBF grad_curr;
        FatropVecBF grad_next;
        FatropVecBF du_inf_curr;
        FatropVecBF du_inf_curr_s;
        // vectors neccessary for inequality constraints
        FatropVecBF s_curr;
        FatropVecBF s_next;
        FatropVecBF delta_s;
        FatropVecBF zL_curr;
        FatropVecBF zL_next;
        FatropVecBF zU_curr;
        FatropVecBF zU_next;
        FatropVecBF delta_zL;
        FatropVecBF delta_zU;
        FatropVecBF s_lower;
        FatropVecBF s_upper;
        struct EvalCache
        {
            struct Instance
            {
                bool evaluated = false;
                double value = 0.0;
                double SetValue(const double value_)
                {
                    value = value_;
                    evaluated = true;
                    return value;
                }
            };
            Instance cv_linf;
            Instance cv_l1;
            Instance lam_linf;
            Instance lam_l1;
            Instance du_inf_linf;
        };
        EvalCache cache_curr;
        EvalCache cache_next;
        double obj_curr;
        double theta_min = 1e-4;
        const RefCountPtr<FatropParams> params;
        // algorithm parameters
        double smax;
        double kappa1;
        double kappa2;
        double kappa_d;
    };
}
#endif // FATROPDATAINCLUDED