/**
 * @file OCPKKT.hpp
 * @author your name (you@domain.com)
 * @brief this file contains classes to represent a KKT matrix by blasfeo submatrices
 * @version 0.1
 * @date 2021-11-10
 *
 * @copyright Copyright (c) 2021
 *
 */
#ifndef FATROP_OCP_KKT_INCLUDED
#define FATROP_OCP_KKT_INCLUDED
#include "../BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "OCPDims.hpp"
#include "../AUX/Aux.hpp"
#include "../AUX/FatropVector.hpp"
using namespace std;
namespace fatrop
{
    /** \brief this class contains all information to represent the KKT system of an equality constrained OCP*/
    class OCPKKTMemory
    {
    public:
        OCPKKTMemory(const OCPDims &dims, MemoryAllocator &fma) : K(dims.K), nu(dims.K, dims.nu, fma), nx(dims.K, dims.nx, fma), ng(dims.K, dims.ng, fma), RSQrqt(dims.nu + dims.nx + 1, dims.nu + dims.nx, dims.K), BAbt(dims.nu + dims.nx + 1, rotate(dims.nx, 1), dims.K), Ggt(dims.nu + dims.nx + 1, dims.ng, dims.K), aux(dims, fma){};

        OCPKKTMemory(const OCPKKTMemory &cpy) = delete;
        int K;
        FatropVector<int> nu;
        FatropVector<int> nx;
        FatropVector<int> ng;
        /// small-scale Hessian
        FatropMemoryMatBF RSQrqt;
        /// small-scale Jacobian dynamics
        FatropMemoryMatBF BAbt;
        /// small-scale Jacobian stagewise eq constraints
        FatropMemoryMatBF Ggt;
        class OCPAux
        {
        public:
            OCPAux(const OCPDims &dims, MemoryAllocator &fma) : ux_offs(dims.K, offsets(dims.nx + dims.nu), fma), g_offs(dims.K, offsets(dims.ng), fma), dyn_eq_offs(dims.K - 1, offsets(rotate(dims.nx, 1)) + sum(dims.ng), fma), max_nu(max(dims.nu)), max_nx(max(dims.nx)), max_ng(max(dims.ng)){};
            /// offset arrays are used for efficiency
            const FatropVector<int> ux_offs;
            /// offset arrays are used for efficiency
            const FatropVector<int> g_offs;
            const FatropVector<int> dyn_eq_offs;
            int max_nu;
            int max_nx;
            int max_ng;
        };
        OCPAux aux;
    };
  

} // namespace fatrop
#endif // FATROP_OCP_KKT_INCLUDED