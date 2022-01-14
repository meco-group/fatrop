// Basic OCP template: initial and terminal constraints eq constraints, Function evaluation provided by Casadi CodeGen
#ifndef OCPTEMPLATEBASICINCLUDED
#define OCPTEMPLATEBASICINCLUDED
#include "OCPTemplate.hpp"
#include <string>
using namespace std;
namespace fatrop
{
    class OCPTemplateBasic: public OCPTemplate
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
        }
        int get_nxk(const int k) const
        {
            return nx_;
        }
        int get_nuk(const int k) const
        {
            return nu_;
        }
        int get_ngk(const int k) const
        {
            if (k == 0)
                return ngI_;
            if (k == K_ - 1)
                return ngF_;
            return 0;
        }
        int get_horizon_length() const { return K_; };
        int eval_BAbtk(const double *states_kp1,
                       const double *scales_states_kp1,
                       const double *states_k,
                       const double *scales_states_k,
                       const double *inputs_k,
                       const double *scales_inputs_k,
                       const double *scales_lam,
                       MAT *res,
                       const int k)
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
                         const int k)
        {
            return 0;
        };
        int eval_Ggtk(const double *states_k,
                      const double *scales_states_k,
                      const double *inputs_k,
                      const double *scales_inputs_k,
                      const double *scales,
                      MAT *res,
                      const int k)
        {
            return 0;
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
