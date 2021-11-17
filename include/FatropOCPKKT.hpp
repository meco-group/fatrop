/**
 * @file FatropOCPKKT.hpp
 * @author your name (you@domain.com)
 * @brief this file contains classes to represent a KKT matrix by blasfeo submatrices 
 * @version 0.1
 * @date 2021-11-10
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef FATROPOCPKKTINDLUCED
#define FATROPOCPKKTINDLUCED
#include "FatropLinearAlgebraBlasfeo.hpp"
#include "FatropOCP.hpp"
#include "FatropAux.hpp"
using namespace std;
namespace fatrop
{
    /** \brief this class contains all information to represent the KKT system of an equality constrained OCP*/
    class OCP_KKT
    {
    public:
        OCP_KKT(const OCP_dims &dims, fatrop_memory_allocator &fma) : nu(dims.K, vector<int>(dims.nu), fma), nx(dims.K, vector<int>(dims.nx), fma), ng(dims.K, vector<int>(dims.ng), fma), RSQrqt(vector<int>(dims.nu + dims.nx + 1), vector<int>(dims.nu + dims.nx), dims.K, fma), BAbt(vector<int>(dims.nu + dims.nx + 1), vector<int>(dims.nx), dims.K, fma), Ggt(vector<int>(dims.nu + dims.nx + 1), vector<int>(dims.ng), dims.K, fma){};
        fatrop_memory_el<int> nu;
        fatrop_memory_el<int> nx;
        fatrop_memory_el<int> ng;
        /// small-scale Hessian
        fatrop_memory_matrix_bf RSQrqt;
        /// small-scale Jacobian dynamics
        fatrop_memory_matrix_bf BAbt;
        /// small-scale Jacobian stagewise eq constraints
        fatrop_memory_matrix_bf Ggt;
    };
    class OCP_aux
    {
    public:
        OCP_aux(const OCP_dims &dims, fatrop_memory_allocator &fma) : ux_offs(dims.K, csum(dims.nx + dims.nu), fma), g_offs(dims.K, csum(dims.ng), fma), max_nu(max(dims.nu)), max_nx(max(dims.nx)), max_ng(max(dims.ng)){};
        /// offset arrays are used for efficiency
        fatrop_memory_el<int> ux_offs;
        /// offset arrays are used for efficiency
        fatrop_memory_el<int> g_offs;
        int max_nu;
        int max_nx;
        int max_ng;
    };
    class OCP_KKT_solver
    {
    public:
        OCP_KKT_solver(const OCP_dims &dims, fatrop_memory_allocator &fma) : aux(dims, fma){};
        /// solve a KKT system
        // void fact_solve(OCP_KKT *KKT, VEC *ux, VEC *lambda)
        // {
        // }
        OCP_aux aux;
    };
} // namespace fatrop
#endif // FATROPOCPKKTINDLUCED