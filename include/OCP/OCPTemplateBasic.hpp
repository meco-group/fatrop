// Basic OCP template: initial and terminal constraints eq constraints, Function evaluation provided by Casadi CodeGen
#ifndef OCPTEMPLATEBASICINCLUDED
#define OCPTEMPLATEBASICINCLUDED
#include "OCPTemplate.hpp"
#include <string>
#include <iostream>
#include <AUX/DynamicLib.hpp>
#include <AUX/SmartPtr.hpp>
#include "FUNCTION_EVALUATION/CasadiCodegen.hpp"
using namespace std;
namespace fatrop
{
    class OCPTemplateBasic : public OCPTemplate
    {
    public:
        OCPTemplateBasic(const int nu,
                         const int nx,
                         const int ngI,
                         const int ngF,
                         const int K,
                         const EvalCasGen &BAbtf,
                         const EvalCasGen &RSQrqtf,
                         const EvalCasGen &RSQrqtFf,
                         const EvalCasGen &GgtIf,
                         const EvalCasGen &GgtFf) : nu_(nu),
                                                    nx_(nx),
                                                    ngI_(ngI),
                                                    ngF_(ngF),
                                                    K_(K),
                                                    BAbtf(BAbtf),
                                                    RSQrqtf(RSQrqtf),
                                                    RSQrqtFf(RSQrqtFf),
                                                    GgtIf(GgtIf),
                                                    GgtFf(GgtFf){};
        static OCPTemplateBasic from_shared_lib(const string &filename, const int K)
        {
            RefCountPtr<DLLoader> handle = new DLLoader(filename);
            EvalCasGen BAbtf(handle, "BAbt");
            EvalCasGen RSQrqtf(handle, "RSQrqt");
            EvalCasGen RSQrqtFf(handle, "RSQrqtF");
            EvalCasGen GgtIf(handle, "GgtI");
            EvalCasGen GgtFf(handle, "GgtF");
            const int nx = BAbtf.out_n;
            const int nu = BAbtf.out_m - nx - 1;
            const int ngI = GgtIf.out_n;
            const int ngF = GgtFf.out_n;
            return OCPTemplateBasic(nu, nx, ngI, ngF, K, BAbtf, RSQrqtf, RSQrqtFf, GgtIf, GgtFf);
        }
        int get_nxk(const int k) const override
        {
            return nx_;
        }
        int get_nuk(const int k) const override
        {
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
            if (k != 0 && k != K_ - 1)
                return 0;
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

    private:
        const int nu_;
        const int nx_;
        const int ngI_;
        const int ngF_;
        const int K_;
        EvalCasGen BAbtf;
        EvalCasGen RSQrqtf;
        EvalCasGen RSQrqtFf;
        EvalCasGen GgtIf;
        EvalCasGen GgtFf;
    };
}
#endif // OCPTEMPLATEBASICINCLUDED
