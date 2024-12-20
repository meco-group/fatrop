#ifndef __fatrop_context_generic_hpp__
#define __fatrop_context_generic_hpp__

// elementary types
namespace fatrop
{
        typedef double Scalar;
        typedef int Index;
} // namespace fatrop

// blasfeo types and functions
#define VEC blasfeo_dvec
#define VECEL BLASFEO_DVECEL
#define AXPBY blasfeo_daxpby
#define AXPY blasfeo_daxpy
#define VECSE blasfeo_dvecse
#define VECSC blasfeo_dvecsc
#define VECCPSC blasfeo_dveccpsc
#define VECCP blasfeo_dveccp
#define VECMUL blasfeo_dvecmul

#define ALLOCATE_VEC blasfeo_allocate_dvec
#define FREE_VEC blasfeo_free_dvec

#endif // __fatrop_context_generic_hpp__