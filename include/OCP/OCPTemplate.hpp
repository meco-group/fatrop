#ifndef OCPTEMPLATEINCLUDED
#define OCPTEMPLATEINCLUDED
#include "OCP/OCPKKT.hpp"
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "AUX/SmartPtr.hpp"
#include "AUX/FatropVector.hpp"
#include "OCPDims.hpp"
namespace fatrop
{
    class OCPTemplate : public RefCountedObj
    {
    public:
        OCPTemplate() : nxexpr(RefCountPtr<OCPTemplate>(this)), nuexpr(RefCountPtr<OCPTemplate>(this)), ngexpr(RefCountPtr<OCPTemplate>(this)){};
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

    private:
        class nxExpr : public VecExpr<nxExpr, int>
        {
        public:
            nxExpr(const RefCountPtr<OCPTemplate> &parent) : parent(parent){};
            int getEl(const int ai) const { return parent->get_nxk(ai); };
            int size() const { return parent->get_horizon_length(); };

        private:
            const RefCountPtr<OCPTemplate> parent;
        };
        class nuExpr : public VecExpr<nxExpr, int>
        {
        public:
            nuExpr(const RefCountPtr<OCPTemplate> &parent) : parent(parent){};
            int getEl(const int ai) const { return parent->get_nuk(ai); };
            int size() const { return parent->get_horizon_length(); };

        private:
            const RefCountPtr<OCPTemplate> parent;
        };
        class ngExpr : public VecExpr<nxExpr, int>
        {
        public:
            ngExpr(const RefCountPtr<OCPTemplate> &parent) : parent(parent){};
            int getEl(const int ai) const { return parent->get_ngk(ai); };
            int size() const { return parent->get_horizon_length(); };

        private:
            const RefCountPtr<OCPTemplate> parent;
        };

    public:
        nxExpr nxexpr;
        nxExpr nuexpr;
        nxExpr ngexpr;
        OCPDims GetDims()
        {
            OCPDims res;
            res.K = get_horizon_length();
            res.nu = nuexpr;
            res.nx = nxexpr;
            res.ng = ngexpr;
            return res;
        }
    };
};     // namespace fatrop
#endif // OCPTEMPLATEINCLUDED