#ifndef __fatrop_context_generic_hpp__
#define __fatrop_context_generic_hpp__

#include "fatrop/common/exception.hpp"
#include "fatrop/linear_algebra/blasfeo_wrapper.hpp"

// elementary types
namespace fatrop
{
        typedef double Scalar;
        typedef int Index;
} // namespace fatrop

// blasfeo types and functions
#define VEC blasfeo_dvec
#define VECEL BLASFEO_DVECEL

// Vector-related BLASFEO definitions
#define AXPBY fatrop::blasfeo_daxpby_debug
#define AXPY fatrop::blasfeo_daxpy_debug
#define VECSE fatrop::blasfeo_dvecse_debug
#define VECSC fatrop::blasfeo_dvecsc_debug
#define VECCPSC fatrop::blasfeo_dveccpsc_debug
#define VECCP fatrop::blasfeo_dveccp_debug
#define VECMUL fatrop::blasfeo_dvecmul_debug

#define ALLOCATE_VEC blasfeo_allocate_dvec
#define FREE_VEC blasfeo_free_dvec

// Matrix-related BLASFEO definitions
#define MAT blasfeo_dmat
#define MATEL BLASFEO_DMATEL
#define ALLOCATE_MAT blasfeo_allocate_dmat
#define FREE_MAT blasfeo_free_dmat

#endif // __fatrop_context_generic_hpp__
