#ifndef OCPTEMPLATEINCLUDED
#define OCPTEMPLATEINCLUDED
#include "OCP/OCPKKT.hpp"
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
namespace fatrop
{
    template<typename Derived>
    class OCPTemplate
    {
    public:
        virtual int get_nxk(const int k);
        virtual int get_nuk(const int k);
        virtual int get_ngk(const int k);
        virtual int get_horizon_lenght();
        virtual int eval_BAbtk(const double ** args, MAT *res, const int k) const = 0;
        virtual int eval_RSQrqtk(const double **args, MAT *res, const int k) const = 0;
        virtual int eval_Ggtk(const double **args, MAT *res, const int k) const = 0;
    };
};     // namespace fatrop
#endif // OCPTEMPLATEINCLUDED