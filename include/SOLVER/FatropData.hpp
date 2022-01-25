// solver data
#ifndef FATROPDATAINCLUDED
#define FATROPDATAINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "TEMPLATES/NLPAlg.hpp"
#include "AUX/Common.hpp"
using namespace std;
namespace fatrop
{
    struct FatropData : public RefCountedObj
    {
        FatropData(const NLPDims &nlpdims) : nlpdims(nlpdims),
                                             memvars(nlpdims.nvars, 6),
                                             memeqs(nlpdims.neqs, 6),
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
                                             grad_next(memvars[5])

        {
        }
        int TryStep(double alpha_primal, double alpha_dual)
        {
            axpy(alpha_primal, delta_x, x_curr, x_next);
            axpby(alpha_primal, lam_calc, 1 - alpha_primal, lam_curr, lam_next);
            // reset evaluation flags
            evalcachenext = EvalCache();
            return 0;
        }
        int TakeStep()
        {
            x_curr.SwapWith(x_next);
            lam_curr.SwapWith(lam_next);
            evalcachecurr = evalcachenext;
            return 0;
        }
        double CVLinfCurr()
        {
            if(evalcachecurr.CVLinfEvaluated) return evalcachecurr.CVLinf;
            evalcachecurr.CVLinfEvaluated = true;
            return Linf(g_curr);
        }
        double CVLinfNext()
        {
            if(evalcachenext.CVLinfEvaluated) return evalcachenext.CVLinf;
            evalcachenext.CVLinfEvaluated = true;
            return Linf(g_next);
        }
        double CVL1Curr()
        {
            if(evalcachecurr.CVL1Evaluated) return evalcachecurr.CVL1;
            evalcachecurr.CVL1Evaluated = true;
            return L1(g_curr);
        }
        double CVL1Next()
        {
            if(evalcachenext.CVL1Evaluated) return evalcachenext.CVL1;
            evalcachenext.CVL1Evaluated = true;
            return L1(g_next);
        }
        const NLPDims nlpdims;
        double obj_scale = 1.0;
        FatropMemoryVecBF memvars;
        FatropMemoryVecBF memeqs;
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
        struct EvalCache
        {
            bool CVLinfEvaluated = false;
            double CVLinf = 0.0;
            bool CVL1Evaluated = false;
            double CVL1 = 0.0;
        };
        EvalCache evalcachecurr; 
        EvalCache evalcachenext; 
    };
}
#endif // FATROPDATAINCLUDED