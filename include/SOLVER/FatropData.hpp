// solver data
#ifndef FATROPDATAINCLUDED
#define FATROPDATAINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "TEMPLATES/NLPAlg.hpp"
#include "AUX/Common.hpp"
#include "FatropParams.hpp"
using namespace std;
namespace fatrop
{
#define CACHEMACRO(instance, val) instance.evaluated ? instance.value : instance.SetValue(val)
    struct FatropData : public RefCountedObj
    {
        FatropData(const NLPDims &nlpdims, const RefCountPtr<FatropParams> &params) : nlpdims(nlpdims),
                                                                                      memvars(nlpdims.nvars, 7),
                                                                                      memeqs(nlpdims.neqs, 6),
                                                                                      memineqs(nlpdims.nineqs, 7),
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
                                                                                      s_curr(memineqs[0]),
                                                                                      s_next(memineqs[1]),
                                                                                      delta_s(memineqs[2]),
                                                                                      zU_curr(memineqs[3]),
                                                                                      zU_next(memineqs[4]),
                                                                                    //   z_next(memineqs[4]),
                                                                                      delta_z(memineqs[5]),
                                                                                      s_upper(memineqs[6]),
                                                                                    //   s_lower(memineqs[6]),
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
            return 0;
        }
        int TryStep(double alpha_primal, double alpha_dual)
        {
            axpy(alpha_primal, delta_x, x_curr, x_next);
            axpby(alpha_primal, lam_calc, 1.0 - alpha_primal, lam_curr, lam_next);
            // reset evaluation flags
            cache_next = EvalCache();
            return 0;
        }
        int TakeStep()
        {
            // TODO make a struct which containts vectors associated with curr <-> next
            x_curr.SwapWith(x_next);
            lam_curr.SwapWith(lam_next);
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

        const NLPDims nlpdims;
        double obj_scale = 1.0;
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
        // vectors neccessary for inequality constraints
        FatropVecBF s_curr;
        FatropVecBF s_next;
        FatropVecBF delta_s;
        FatropVecBF zU_curr;
        FatropVecBF zU_next;
        FatropVecBF delta_z;
        FatropVecBF s_upper;
        // FatropVecBF s_upper;
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