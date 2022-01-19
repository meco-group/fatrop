#ifndef OCPTEMPLATEINCLUDED
#define OCPTEMPLATEINCLUDED
#include "OCP/OCPKKT.hpp"
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "AUX/SmartPtr.hpp"
#include "AUX/FatropVector.hpp"
#include "OCPDims.hpp"
namespace fatrop
{
    class BFOCP : public RefCountedObj
    {
    public:
        virtual int get_nxk(const int k) const = 0;
        virtual int get_nuk(const int k) const = 0;
        virtual int get_ngk(const int k) const = 0;
        virtual int get_horizon_length() const = 0;
        virtual int eval_BAbtk(const double *states_kp1,
                               const double *scales_states_kp1,
                               const double *states_k,
                               const double *scales_states_k,
                               const double *inputs_k,
                               const double *scales_inputs_k,
                               const double *scales_lam,
                               MAT *res,
                               const int k) = 0;
        virtual int eval_RSQrqtk(const double *objective_scale,
                                 const double *states_k,
                                 const double *scales_states_k,
                                 const double *inputs_k,
                                 const double *scales_inputs_k,
                                 const double *lam_dyn_k,
                                 const double *scales_lam_dyn_k,
                                 const double *lam_eq_k,
                                 const double *scales_lam_eq_k,
                                 MAT *res,
                                 const int k) = 0;
        virtual int eval_Ggtk(const double *states_k,
                              const double *scales_states_k,
                              const double *inputs_k,
                              const double *scales_inputs_k,
                              const double *scales,
                              MAT *res,
                              const int k) = 0;

    };
};     // namespace fatrop
#endif // OCPTEMPLATEINCLUDED