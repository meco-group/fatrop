#ifndef OCPINCLUDED
#define OCPINCLUDED
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "OCPDims.hpp"
namespace fatrop
{
    class OCP
    {
        public:
        virtual int evalHess(
            OCPKKTMemory* OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &scales_primal_vars,
            const FatropVecBF &lam,
            const FatropVecBF &scales_lam) = 0;
        virtual int evalJac(
            OCPKKTMemory* OCP,
            const FatropVecBF &primal_vars,
            const FatropVecBF &scales_primal_vars,
            const FatropVecBF &scales_lam) = 0;
        virtual OCPDims GetOCPDims() const  = 0;
    };
} // namespace fatrop
#endif // OCPINCLUDED