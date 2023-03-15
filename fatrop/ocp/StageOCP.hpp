// Basic OCP template: initial and terminal constraints eq constraints, Function evaluation provided by Casadi CodeGen
#ifndef OCPTEMPLATEBASICINCLUDED
#define OCPTEMPLATEBASICINCLUDED
#include "OCPAbstract.hpp"
#include <string>
#include <iostream>
#include <aux/DynamicLib.hpp>
#include <aux/SmartPtr.hpp>
#include "function_evaluation/CasadiCodegen.hpp"
#include "json/json.h"
#include <fstream>
#include <sstream>

using namespace std;
namespace fatrop
{
    /* TODO: at this point the StageOCP's implementation implements function evaluation through Casadi Codegen 
    *  it should become an interface that can be implemented by user defined function evaluation methods 
    *  the StageOCP is then characterised by its structure that is a bit less general than the OCPAbstract interface
    *  more in particular it seperates initial, intermediate and terminal stages, where it has different constraints and objectives
    * 
    *  the casadi implementation should become a derived class of this interface
    */ 

    /// BasicOCP is a OCPAbstract implementation that is built on top of CasAdi Codegenerated functions 
    /// it seperates the initial and terminal stages and the intermediate stages, each of these can have different constraints and running objective
    /// the dynamics is the same for all stages
    /// this problem structure corresponds to a single-stage rockit problem, which it is intented to be used for
    class StageOCP : public OCPAbstract
    {
    public:
        StageOCP(const int nu,
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

    public:
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

    private:
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

    class StageOCPBuilder
    {
    public:
        static shared_ptr<StageOCP> FromRockitInterface(const shared_ptr<DLHandler> &handle, const json::jobject& json_spec)
        {

            // set up ocp
            const bool GN = false;
            // shared_ptr<DLHandler> handle = make_shared<DLHandler>(functions);
            // std::ifstream t(json_spec_file);
            // std::stringstream buffer;
            // buffer << t.rdbuf();
            // json::jobject json_spec = json::jobject::parse(buffer.str());
            const int K = json_spec["K"];
            const int nx = json_spec["nx"];
            const int nu = json_spec["nu"];
            const int ngI = json_spec["ngI"];
            const int ng = json_spec["ng"];
            const int ngF = json_spec["ngF"];
            const int ng_ineqI = json_spec["ng_ineqI"];
            const int ng_ineq = json_spec["ng_ineq"];
            const int ng_ineqF = json_spec["ng_ineqF"];
            const int no_stage_params = json_spec["n_stage_params"];
            const int no_global_params = json_spec["n_global_params"];
            vector<double> lowerI = json_spec["lowerI"].get_number_array<double>("%lf");
            vector<double> upperI = json_spec["upperI"].get_number_array<double>("%lf");
            vector<double> lower = json_spec["lower"].get_number_array<double>("%lf");
            vector<double> upper = json_spec["upper"].get_number_array<double>("%lf");
            vector<double> lowerF = json_spec["lowerF"].get_number_array<double>("%lf");
            vector<double> upperF = json_spec["upperF"].get_number_array<double>("%lf");
            lower.insert(lower.begin(), lowerI.begin(), lowerI.end());
            upper.insert(upper.begin(), upperI.begin(), upperI.end());
            lower.insert(lower.end(), lowerF.begin(), lowerF.end());
            upper.insert(upper.end(), upperF.begin(), upperF.end());
            vector<double> initial_u = json_spec["initial_u"].get_number_array<double>("%lf");
            vector<double> initial_x = json_spec["initial_x"].get_number_array<double>("%lf");
            EvalCasGen BAbtf(handle, "BAbt");
            EvalCasGen bkf(handle, "bk");
            EvalCasGen RSQrqtIf = GN ? EvalCasGen(handle, "RSQrqtIGN") : EvalCasGen(handle, "RSQrqtI");
            EvalCasGen rqIf(handle, "rqI");
            EvalCasGen RSQrqtf = GN ? EvalCasGen(handle, "RSQrqtGN") : EvalCasGen(handle, "RSQrqt");
            EvalCasGen rqf(handle, "rqk");
            EvalCasGen RSQrqtFf = GN ? EvalCasGen(handle, "RSQrqtFGN") : EvalCasGen(handle, "RSQrqtF");
            EvalCasGen rqFf(handle, "rqF");
            EvalCasGen GgtIf(handle, "GgtI");
            EvalCasGen gIf(handle, "gI");
            EvalCasGen Ggtf(handle, "Ggt");
            EvalCasGen gf(handle, "g");
            EvalCasGen GgtFf(handle, "GgtF");
            EvalCasGen gFf(handle, "gF");
            EvalCasGen LIf(handle, "LI");
            EvalCasGen Lkf(handle, "Lk");
            EvalCasGen LFf(handle, "LF");
            EvalCasGen GgineqItf(handle, "GgineqIt");
            EvalCasGen gineqIf(handle, "gineqI");
            EvalCasGen Ggineqtf(handle, "Ggineqt");
            EvalCasGen gineqf(handle, "gineq");
            EvalCasGen GgineqFtf(handle, "GgineqFt");
            EvalCasGen gineqFf(handle, "gineqF");
            shared_ptr<StageOCP> stageocp = make_shared<StageOCP>(nu, nx, ngI, ng, ngF, ng_ineqI, ng_ineq, ng_ineqF, no_stage_params, no_global_params, K,
                                                                  BAbtf,
                                                                  bkf,
                                                                  RSQrqtIf,
                                                                  rqIf,
                                                                  RSQrqtf,
                                                                  rqf,
                                                                  RSQrqtFf,
                                                                  rqFf,
                                                                  GgtIf,
                                                                  gIf,
                                                                  Ggtf,
                                                                  gf,
                                                                  GgtFf,
                                                                  gFf,
                                                                  GgineqItf,
                                                                  gineqIf,
                                                                  Ggineqtf,
                                                                  gineqf,
                                                                  GgineqFtf,
                                                                  gineqFf,
                                                                  LIf,
                                                                  Lkf,
                                                                  LFf, lower, upper,
                                                                  json_spec["stage_params"].get_number_array<double>("%lf"),
                                                                  json_spec["global_params"].get_number_array<double>("%lf"),
                                                                  initial_u, initial_x);
            return stageocp;
        }
    };
}
#endif // OCPTEMPLATEBASICINCLUDED