#ifndef __fatrop_context_generic_hpp__
#define __fatrop_context_generic_hpp__

#include "fatrop/common/exception.hpp"

// elementary types
namespace fatrop
{
    typedef double Scalar;
    typedef int Index;
} // namespace fatrop


// blasfeo types and functions
#define VEC blasfeo_dvec
#define VECEL BLASFEO_DVECEL
#define MAT blasfeo_dmat
#define MATEL BLASFEO_DMATEL


// Vector-related BLASFEO definitions
#define AXPBY blasfeo_daxpby
#define AXPY blasfeo_daxpy
#define VECPE blasfeo_dvecpe
#define VECPEI blasfeo_dvecpei
#define VECSE blasfeo_dvecse
#define VECSC blasfeo_dvecsc
#define VECCPSC blasfeo_dveccpsc
#define VECCP blasfeo_dveccp
#define VECMUL blasfeo_dvecmul
#define DOT blasfeo_ddot

#define ALLOCATE_VEC blasfeo_allocate_dvec
#define FREE_VEC blasfeo_free_dvec
#define MEMSIZE_VEC blasfeo_memsize_dvec

// Matrix-related BLASFEO definitions
#define ROWPE blasfeo_drowpe
#define ROWPEI blasfeo_drowpei
#define COLPE blasfeo_dcolpe
#define COLPEI blasfeo_dcolpei
#define ROWSW blasfeo_drowsw
#define COLSW blasfeo_dcolsw
#define GEAD blasfeo_dgead
#define GECP blasfeo_dgecp
#define GESC blasfeo_dgesc
#define TRSM_RLTN blasfeo_dtrsm_rltn
#define TRSM_RLNN blasfeo_dtrsm_rlnn
#define GEMM_NT blasfeo_dgemm_nt
#define GEAD blasfeo_dgead
#define SYRK_LN_MN blasfeo_dsyrk_ln_mn
#define SYRK_LN blasfeo_dsyrk_ln
#define GETR blasfeo_dgetr
#define TRTR_L blasfeo_dtrtr_l
#define POTRF_L_MN blasfeo_dpotrf_l_mn
#define ROWEX blasfeo_drowex
#define ROWIN blasfeo_drowin
#define COLIN blasfeo_dcolin
#define TRSV_LTN blasfeo_dtrsv_ltn
#define TRSV_LNN blasfeo_dtrsv_lnn
#define TRSV_UTN blasfeo_dtrsv_utn
#define GEMV_T blasfeo_dgemv_t
#define GEMV_N blasfeo_dgemv_n
#define GESE blasfeo_dgese
#define DIARE blasfeo_ddiare
#define DIAAD blasfeo_ddiaad
#define COLSC blasfeo_dcolsc
#define VECMULACC blasfeo_dvecmulacc
#define GER blasfeo_dger

#define PACK_MAT blasfeo_pack_dmat
#define UNPACK_MAT blasfeo_unpack_dmat
#define PACK_VEC blasfeo_pack_dvec
#define UNPACK_VEC blasfeo_unpack_dvec

#define ALLOCATE_MAT blasfeo_allocate_dmat
#define FREE_MAT blasfeo_free_dmat
#define MEMSIZE_MAT blasfeo_memsize_dmat

#endif // __fatrop_context_generic_hpp__
