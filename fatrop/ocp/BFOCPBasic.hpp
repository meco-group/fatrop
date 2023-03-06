// Basic OCP template: initial and terminal constraints eq constraints, Function evaluation provided by Casadi CodeGen
#ifndef OCPTEMPLATEBASICINCLUDED
#define OCPTEMPLATEBASICINCLUDED
#include "OCPAbstract.hpp"
#include <string>
#include <iostream>
#include <aux/DynamicLib.hpp>
#include <aux/SmartPtr.hpp>
#include "function_evaluation/CasadiCodegen.hpp"
using namespace std;
namespace fatrop
{
    class BFOCPBasic : public OCPAbstract
    {
    public:
        BFOCPBasic(const int nu,
                   const int nx,
                   const int ngI,
                   const int ng,
                   const int ngF,
                   const int ng_ineqI,
                   const int ng_ineq,
                   const int ng_ineqF,
                   const int n_stage_params,
                   const int n_global_params,
                   const int K,
                   const EvalCasGen &BAbtf,
                   const EvalCasGen &bkf,
                   const EvalCasGen &RSQrqtIf,
                   const EvalCasGen &rqIf,
                   const EvalCasGen &RSQrqtf,
                   const EvalCasGen &rqf,
                   const EvalCasGen &RSQrqtFf,
                   const EvalCasGen &rqFf,
                   const EvalCasGen &GgtIf,
                   const EvalCasGen &gIf,
                   const EvalCasGen &Ggtf,
                   const EvalCasGen &gf,
                   const EvalCasGen &GgtFf,
                   const EvalCasGen &gFf,
                   const EvalCasGen &Ggt_ineqIf,
                   const EvalCasGen &gineqIf,
                   const EvalCasGen &Ggt_ineqf,
                   const EvalCasGen &gineqf,
                   const EvalCasGen &Ggt_ineqFf,
                   const EvalCasGen &gineqFf,
                   const EvalCasGen &LkIf,
                   const EvalCasGen &Lkf,
                   const EvalCasGen &LFf,
                   const vector<double> &bounds_L,
                   const vector<double> &bounds_U, 
                   const vector<double> &stage_params,
                   const vector<double> &global_params,
                   const vector<double> &initial_u,
                   const vector<double> &initial_x);
        int get_nxk(const int k) const override;
        int get_nuk(const int k) const override;
        int get_ngk(const int k) const override;
        int get_ng_ineq_k(const int k) const override;
        int get_n_global_params() const;
        int get_n_stage_params_k(const int k) const override;
        int get_default_stage_paramsk(double *stage_params_res, const int k) const override
        {
            int offs = k * n_stage_params_;
            const double *stage_params_p = stage_params.data();
            for (int i = 0; i < n_stage_params_; i++)
            {
                stage_params_res[i] = stage_params_p[offs + i];
            }
            return 0;
        }
        int get_default_global_params(double *global_params_res) const override
        {
            const double *global_params_p = global_params.data();
            for (int i = 0; i < n_global_params_; i++)
            {
                global_params_res[i] = global_params_p[i];
            }
            return 0;
        }
        int get_horizon_length() const;
        int eval_BAbtk(const double *states_kp1,
                       const double *inputs_k,
                       const double *states_k,
                       const double *stage_params_k,
                       const double *global_params,
                       MAT *res,
                       const int k) override;
        int eval_RSQrqtk(const double *objective_scale,
                         const double *inputs_k,
                         const double *states_k,
                         const double *lam_dyn_k,
                         const double *lam_eq_k,
                         const double *lam_ineq_k,
                         const double *stage_params_k,
                         const double *global_params,
                         MAT *res,
                         const int k) override;
        int eval_Ggtk(
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            MAT *res,
            const int k) override;
        int eval_Ggt_ineqk(
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            MAT *res,
            const int k) override;
        int eval_bk(
            const double *states_kp1,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            double *constraint_violation_k,
            const int k) override;
        int eval_gk(
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            double *res,
            const int k) override;
        int eval_gineqk(
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            double *res,
            const int k) override;
        int eval_rqk(
            const double *objective_scale,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            double *res,
            const int k) override;

        int eval_Lk(
            const double *objective_scale,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            double *res,
            const int k) override;
        int get_initial_xk(double *xk, const int k) const override
        {
            const double *initial_x_p = initial_x.data();
            for (int i = 0; i < nx_; i++)
            {
                xk[i] = initial_x_p[i + k * nx_];
            }
            return 0;
        };
        int get_initial_uk(double *uk, const int k) const override
        {
            const double *initial_u_p = initial_u.data();
            for (int i = 0; i < nu_; i++)
            {
                uk[i] = initial_u_p[i + k * nu_];
            }
            return 0;
        };
        int set_initial_xk(double *xk, const int k)
        {
            double *initial_x_p = initial_x.data();
            for (int i = 0; i < nx_; i++)
            {
                initial_x_p[i + k * nx_] = xk[i];
            }
            return 0;
        };
        int set_initial_uk(double *uk, const int k)
        {
            double *initial_u_p = initial_u.data();
            for (int i = 0; i < nu_; i++)
            {
                initial_u_p[i + k * nu_] = uk[i];
            }
            return 0;
        };
        int get_boundsk(double *lower, double *upper, const int k) const override
        {
            const double *bounds_L_p = bounds_L.data();
            const double *bounds_U_p = bounds_U.data();
            int offs = 0;
            int ngineq = ng_ineq_;
            if (k == 0)
            {
                offs = 0;
            }
            else
            {
                offs = ng_ineqI_ + (k - 1) * ng_ineq_;
            }
            if (k == 0)
            {
                ngineq = ng_ineqI_;
            }
            else if (k == K_ - 1)
            {
                ngineq = ng_ineqF_;
            }

            for (int i = 0; i < ngineq; i++)
            {
                lower[i] = bounds_L_p[i + offs];
                upper[i] = bounds_U_p[i + offs];
            }
            return 0;
        };

    private:
        const int nu_;
        const int nx_;
        const int ngI_;
        const int ng_;
        const int ngF_;
        const int ng_ineqI_;
        const int ng_ineq_;
        const int ng_ineqF_;
        const int n_stage_params_;
        const int n_global_params_;
        const int K_;
        EvalCasGen BAbtf;
        EvalCasGen bkf;
        EvalCasGen RSQrqtIf;
        EvalCasGen rqIf;
        EvalCasGen RSQrqtf;
        EvalCasGen rqf;
        EvalCasGen RSQrqtFf;
        EvalCasGen rqFf;
        EvalCasGen GgtIf;
        EvalCasGen gIf;
        EvalCasGen Ggtf;
        EvalCasGen gf;
        EvalCasGen GgtFf;
        EvalCasGen gFf;
        EvalCasGen Ggt_ineqIf;
        EvalCasGen g_ineqIf;
        EvalCasGen Ggt_ineqf;
        EvalCasGen g_ineqf;
        EvalCasGen Ggt_ineqFf;
        EvalCasGen g_ineqFf;
        EvalCasGen LkIf;
        EvalCasGen Lkf;
        EvalCasGen LFf;
        vector<double> initial_x;
        vector<double> initial_u;
        vector<double> bounds_L;
        vector<double> bounds_U;
        vector<double> stage_params;
        vector<double> global_params;
    };
}
#endif // OCPTEMPLATEBASICINCLUDED