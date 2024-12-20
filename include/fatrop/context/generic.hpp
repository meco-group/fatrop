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
#define AXPBY fatrop::blasfeo_daxpby_debug
#define AXPY fatrop::blasfeo_daxpy_debug
#define VECPE fatrop::blasfeo_dvecpe_debug
#define VECPEI fatrop::blasfeo_dvecpei_debug
#define VECSE fatrop::blasfeo_dvecse_debug
#define VECSC fatrop::blasfeo_dvecsc_debug
#define VECCPSC fatrop::blasfeo_dveccpsc_debug
#define VECCP fatrop::blasfeo_dveccp_debug
#define VECMUL fatrop::blasfeo_dvecmul_debug
#define DOT fatrop::blasfeo_ddot_debug

#define ALLOCATE_VEC blasfeo_allocate_dvec
#define FREE_VEC blasfeo_free_dvec

// Matrix-related BLASFEO definitions
#define ROWPE fatrop::blasfeo_drowpe_debug
#define ROWPEI fatrop::blasfeo_drowpei_debug
#define COLPE fatrop::blasfeo_dcolpe_debug
#define COLPEI fatrop::blasfeo_dcolpei_debug
#define ROWSW fatrop::blasfeo_drowsw_debug
#define COLSW fatrop::blasfeo_dcolsw_debug
#define GEAD fatrop::blasfeo_dgead_debug
#define GECP fatrop::blasfeo_dgecp_debug
#define GESC fatrop::blasfeo_dgesc_debug
#define TRSM_RLTN fatrop::blasfeo_dtrsm_rltn_debug
#define GEMM_NT fatrop::blasfeo_dgemm_nt_debug
#define GEAD fatrop::blasfeo_dgead_debug
#define SYRK_LN_MN fatrop::blasfeo_dsyrk_ln_mn_debug
#define SYRK_LN fatrop::blasfeo_dsyrk_ln_debug
#define GETR fatrop::blasfeo_dgetr_debug
#define TRTR_L fatrop::blasfeo_dtrtr_l_debug
#define POTRF_L_MN fatrop::blasfeo_dpotrf_l_mn_debug
#define ROWEX fatrop::blasfeo_drowex_debug
#define ROWIN fatrop::blasfeo_drowin_debug
#define COLIN fatrop::blasfeo_dcolin_debug
#define TRSV_LTN fatrop::blasfeo_dtrsv_ltn_debug
#define TRSV_LNN fatrop::blasfeo_dtrsv_lnn_debug
#define TRSV_UTN fatrop::blasfeo_dtrsv_utn_debug
#define GEMV_T fatrop::blasfeo_dgemv_t_debug
#define GEMV_N fatrop::blasfeo_dgemv_n_debug
#define GESE fatrop::blasfeo_dgese_debug
#define DIARE fatrop::blasfeo_ddiare_debug
#define COLSC fatrop::blasfeo_dcolsc_debug
#define VECMULACC fatrop::blasfeo_dvecmulacc_debug
#define GER fatrop::blasfeo_dger_debug

#define ALLOCATE_MAT blasfeo_allocate_dmat
#define FREE_MAT blasfeo_free_dmat

#endif // __fatrop_context_generic_hpp__
