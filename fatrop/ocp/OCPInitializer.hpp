#ifndef OCPINITIALIZERINCLUDED
#define OCPINITIALIZERINCLUDED
#include "OCPKKT.hpp"
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP->name)
#define AUXMACRO(type, name, suffix) type name##suffix = ((type)OCP->aux.name)
#define SOLVERMACRO(type, name, suffix) type name##suffix = ((type)name)
namespace fatrop
{
    class OCPInitializer
    {
    public:
        /** \brief this method adapts KKT system for initialization, JAC and GRAD are assumed evaluated !! */
        int AdaptKKTInitial(
            OCPKKTMemory *OCP,
            const FatropVecBF &grad,
            FatropVecBF& s);
    };
} // namespace fatrop
#endif //  OCPINITIALIZERINCLUDED