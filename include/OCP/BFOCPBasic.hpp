// Basic OCP template: initial and terminal constraints eq constraints, Function evaluation provided by Casadi CodeGen
#ifndef OCPTEMPLATEBASICINCLUDED
#define OCPTEMPLATEBASICINCLUDED
#include "BFOCP.hpp"
#include <string>
#include <iostream>
#include <AUX/DynamicLib.hpp>
#include <AUX/SmartPtr.hpp>
#include "FUNCTION_EVALUATION/CasadiCodegen.hpp"
using namespace std;
namespace fatrop
{
    class BFOCPBasic : public BFOCP
    {
    public:
        BFOCPBasic(const int nu,
                   const int nx,
                   const int ngI,
                   const int ngF,
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
                   const EvalCasGen &GgtFf,
                   const EvalCasGen &gFf) : nu_(nu),
                                            nx_(nx),
                                            ngI_(ngI),
                                            ngF_(ngF),
                                            K_(K),
                                            BAbtf(BAbtf),
                                            bkf(bkf),
                                            RSQrqtIf(RSQrqtIf),
                                            rqIf(rqIf),
                                            RSQrqtf(RSQrqtf),
                                            rqf(rqf),
                                            RSQrqtFf(RSQrqtFf),
                                            rqFf(rqFf),
                                            GgtIf(GgtIf),
                                            gIf(gIf),
                                            GgtFf(GgtFf),
                                            gFf(gFf) {}
        // TODO Create Builder class
        static BFOCPBasic from_shared_lib(const string &filename, const int K)
        {
            RefCountPtr<DLHandler> handle = new DLHandler(filename);
            EvalCasGen BAbtf(handle, "BAbt");
            EvalCasGen bkf(handle, "bk");
            EvalCasGen RSQrqtIf(handle, "RSQrqtI");
            EvalCasGen rqIf(handle, "rqI");
            EvalCasGen RSQrqtf(handle, "RSQrqt");
            EvalCasGen rqf(handle, "rqk");
            EvalCasGen RSQrqtFf(handle, "RSQrqtF");
            EvalCasGen rqFf(handle, "rqF");
            EvalCasGen GgtIf(handle, "GgtI");
            EvalCasGen gIf(handle, "gI");
            EvalCasGen GgtFf(handle, "GgtF");
            EvalCasGen gFf(handle, "gF");

            const int nx = BAbtf.out_n;
            const int nu = BAbtf.out_m - nx - 1;
            const int ngI = GgtIf.out_n;
            const int ngF = GgtFf.out_n;
            return BFOCPBasic(nu, nx, ngI, ngF, K,
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
                              GgtFf,
                              gFf);
        }
        int get_nxk(const int k) const override
        {
            return nx_;
        }
        int get_nuk(const int k) const override
        {
            if (k == K_ - 1)
                return 0;
            return nu_;
        }
        int get_ngk(const int k) const override
        {
            if (k == 0)
                return ngI_;
            if (k == K_ - 1)
                return ngF_;
            return 0;
        }
        int get_horizon_length() const override { return K_; };
        int eval_BAbtk(const double *states_kp1,
                       const double *scales_states_kp1,
                       const double *states_k,
                       const double *scales_states_k,
                       const double *inputs_k,
                       const double *scales_inputs_k,
                       const double *scales_lam,
                       MAT *res,
                       const int k) override
        {
            const double *args[7];
            args[0] = states_kp1;
            args[1] = scales_states_kp1;
            args[2] = states_k;
            args[3] = scales_states_k;
            args[4] = inputs_k;
            args[5] = scales_inputs_k;
            args[6] = scales_lam;
            return BAbtf.eval_bf(args, res);
        }
        int eval_RSQrqtk(const double *objective_scale,
                         const double *states_k,
                         const double *scales_states_k,
                         const double *inputs_k,
                         const double *scales_inputs_k,
                         const double *lam_dyn_k,
                         const double *scales_lam_dyn_k,
                         const double *lam_eq_k,
                         const double *scales_lam_eq_k,
                         MAT *res,
                         const int k) override
        {
            const double *args[9];
            args[0] = objective_scale;
            args[1] = states_k;
            args[2] = scales_states_k;
            args[3] = inputs_k;
            args[4] = scales_inputs_k;
            args[5] = lam_dyn_k;
            args[6] = scales_lam_dyn_k;
            args[7] = lam_eq_k;
            args[8] = scales_lam_eq_k;
            if (k == 0)
                return RSQrqtIf.eval_bf(args, res);
            if (k == K_ - 1)
                return RSQrqtFf.eval_bf(args, res);
            return RSQrqtf.eval_bf(args, res);
        };
        int eval_Ggtk(const double *states_k,
                      const double *scales_states_k,
                      const double *inputs_k,
                      const double *scales_inputs_k,
                      const double *scales,
                      MAT *res,
                      const int k) override
        {
            const double *args[5];
            args[0] = states_k;
            args[1] = scales_states_k;
            args[2] = inputs_k;
            args[3] = scales_inputs_k;
            args[4] = scales;
            if (k == K_ - 1)
                return GgtFf.eval_bf(args, res);
            if (k == 0)
                return GgtIf.eval_bf(args, res);
            return 1;
        };
        int eval_bk(
            const double *states_kp1,
            const double *scales_states_kp1,
            const double *states_k,
            const double *scales_states_k,
            const double *inputs_k,
            const double *scales_inputs_k,
            const double *scales_lam,
            double *constraint_violation_k,
            const int k) override
        {
            const double *args[7];
            args[0] = states_kp1;
            args[1] = scales_states_kp1;
            args[2] = states_k;
            args[3] = scales_states_k;
            args[4] = inputs_k;
            args[5] = scales_inputs_k;
            args[6] = scales_lam;
            return bkf.eval_array(args, constraint_violation_k);
        };
        virtual int eval_gk(
            const double *states_k,
            const double *scales_states_k,
            const double *inputs_k,
            const double *scales_inputs_k,
            const double *scales,
            double *res,
            const int k) override
        {
            const double *args[5];
            args[0] = states_k;
            args[1] = scales_states_k;
            args[2] = inputs_k;
            args[3] = scales_inputs_k;
            args[4] = scales;
            if (k == K_ - 1)
                return gFf.eval_array(args, res);
            if (k == 0)
                return gIf.eval_array(args, res);
            return 1;
        }
        int eval_rqk(
            const double *objective_scale,
            const double *states_k,
            const double *scales_states_k,
            const double *inputs_k,
            const double *scales_inputs_k,
            double *res,
            const int k) override
        {
            const double *args[7];
            args[0] = objective_scale;
            args[1] = states_k;
            args[2] = scales_states_k;
            args[3] = inputs_k;
            args[4] = scales_inputs_k;
            if (k == K_ - 1)
                return rqFf.eval_array(args, res);
            if (k == 0)
                return rqIf.eval_array(args, res);
            return rqf.eval_array(args, res);
        };

    private:
        const int nu_;
        const int nx_;
        const int ngI_;
        const int ngF_;
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
        EvalCasGen GgtFf;
        EvalCasGen gFf;
    };
}
#endif // OCPTEMPLATEBASICINCLUDED
