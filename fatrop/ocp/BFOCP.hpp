#ifndef OCPTEMPLATEINCLUDED
#define OCPTEMPLATEINCLUDED
#include "ocp/OCPKKT.hpp"
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "aux/SmartPtr.hpp"
#include "aux/FatropVector.hpp"
#include "OCPDims.hpp"
namespace fatrop
{
    class BFOCP 
    {
    public:
        virtual int get_nxk(const int k) const = 0;
        virtual int get_nuk(const int k) const = 0;
        virtual int get_ngk(const int k) const = 0;
        virtual int get_n_stage_params_k(const int k) const = 0;
        virtual int get_n_global_parmas() const = 0;
        virtual int get_ng_ineq_k(const int k) const = 0;
        virtual int get_horizon_length() const = 0;
        virtual int eval_BAbtk(
            const double *states_kp1,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params_k,
            MAT *res,
            const int k) = 0;
        virtual int eval_RSQrqtk(
            const double *objective_scale,
            const double *inputs_k,
            const double *states_k,
            const double *lam_dyn_k,
            const double *lam_eq_k,
            const double *lam_eq_ineq_k,
            const double *stage_params_k,
            const double *global_params_k,
            MAT *res,
            const int k) = 0;
        virtual int eval_Ggtk(
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params_k,
            MAT *res,
            const int k) = 0;
        virtual int eval_Ggt_ineqk(
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params_k,
            MAT *res,
            const int k) = 0;
        virtual int eval_bk(
            const double *states_kp1,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params_k,
            double *res,
            const int k) = 0;
        virtual int eval_gk(
            const double *states_k,
            const double *inputs_k,
            const double *stage_params_k,
            const double *global_params_k,
            double *res,
            const int k) = 0;
        virtual int eval_gineqk(
            const double *states_k,
            const double *inputs_k,
            const double *stage_params_k,
            const double *global_params_k,
            double *res,
            const int k) = 0;
        virtual int eval_rqk(
            const double *objective_scale,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params_k,
            double *res,
            const int k) = 0;
        virtual int eval_Lk(
            const double *objective_scale,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params_k,
            double * res,
            const int k) = 0;
    };
};     // namespace fatrop
#endif // OCPTEMPLATEINCLUDED