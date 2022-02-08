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
        }
        int Reset()
        {
            cache_curr = EvalCache();
            cache_next = EvalCache();
            return 0;
        }
        double EMuCurr(double mu)
        {
            double res = 0.0;
            double lammean = LamMeanCurr();
            double cv = CVLinfCurr();
            double du = DuInfLinfCurr();
            double sd = 0.0;
            if (lammean > smax)
            {
                sd = lammean / smax;
                du /= sd;
            }
            res = MAX(cv, du);
            return res;
        };
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
                    alpha_max_pr = delta_s_i < 0 ? MIN(alpha_max_pr, -tau * (VECEL(s_curr_p, i)-VECEL(s_lower_p, i)) / delta_s_i) : alpha_max_pr;
                    // dual
                    alpha_max_du = delta_Z_i < 0 ? MIN(alpha_max_du, -tau * (VECEL(zL_curr_p, i)) / delta_Z_i) : alpha_max_du;
                }
                if (!isinf(VECEL(s_upper_p, i)))
                {
                    double delta_s_i = VECEL(delta_s_p, i);
                    double delta_Z_i = VECEL(delta_zU_p, i);
                    // primal
                    alpha_max_pr = delta_s_i < 0 ? MIN(alpha_max_pr, -tau * (VECEL(s_upper_p, i) - VECEL(s_curr_p, i)) / delta_s_i) : alpha_max_pr;
                    // dual
                    alpha_max_du = delta_Z_i < 0 ? MIN(alpha_max_du, -tau * (VECEL(zU_curr_p, i)) / delta_Z_i) : alpha_max_du;
                }
            }
            return;
        }

        const NLPDims nlpdims;
        double obj_scale = 1.0;
        int n_ineqs;
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
    };
}
#endif // FATROPDATAINCLUDED