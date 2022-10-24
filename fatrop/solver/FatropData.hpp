// solver data
#ifndef FATROPDATAINCLUDED
#define FATROPDATAINCLUDED
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "templates/NLPAlg.hpp"
#include "aux/Common.hpp"
#include "FatropParams.hpp"
#include <cmath>
#include <memory>
using namespace std;
namespace fatrop
{
#define CACHEMACRO(instance, val) instance.evaluated ? instance.value : instance.SetValue(val)
    struct FatropData  
    {
        FatropData(const NLPDims &nlpdims, const shared_ptr<FatropParams> &params) ;
        void Initialize();
        int Reset();
        int ResetCaches();
        double EMuCurr(double mu);
        int EvalDuInfSlacksEqs();
        double EvalCompSlackInf(double mu);
        double EvalBarrier(double mu, VEC *s_p);
        double EvalBarrierCurr(double mu);
        double EvalBarrierNext(double mu);
        double EvalBarrierLinDecr(double mu);
        int BoundSlacks();
        int AdaptDualBounds(double mu);
        int AcceptInitialization();
        int TryStep(double alpha_primal, double alpha_dual);
        int TakeStep();
        double CVLinfCurr();
        double CVLinfNext();
        double CVL1Curr();
        double CVL1Next();
        double LamL1Curr();
        double LamLinfCurr();
        double LamMeanCurr();
        double ZL1Curr();
        int nIneqsR();
        double LamLinfCalc();
        double DuInfLinfCurr();
        double LinDecrCurr();
        void AlphaMax(double &alpha_max_pr, double &alpha_max_du, double tau);
        void SetBounds(const vector<double>& lowerin, const vector<double>& upperin);
        void RelaxBounds();
        void RelaxBoundsVar(double mu);

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
        FatropVecBF x_initial;
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
        FatropVecBF residual_vec;
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
        FatropVecBF s_lower_orig;
        FatropVecBF s_upper_orig;
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
        const shared_ptr<FatropParams> params;
        // algorithm parameters
        double smax;
        double kappa1;
        double kappa2;
        double kappa_d;
        double kappa_sigma;
        double bound_relax_factor;
        double constr_viol_tol;
    };
}
#endif // FATROPDATAINCLUDED