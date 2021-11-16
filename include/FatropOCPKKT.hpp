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
    class OCP_KKT
    {
    public:
        OCP_KKT(const fatrop_OCP_dims &dims, fatrop_memory_allocator &fma) : nu(dims.K, vector<int>(dims.nu), fma), nx(dims.K, vector<int>(dims.nx), fma), ng(dims.K,vector<int>(dims.ng), fma), RSQrqt(vector<int>(dims.nu + dims.nx + 1), vector<int>(dims.nu + dims.nx), dims.K, fma), BAbt(vector<int>(dims.nu + dims.nx + 1), vector<int>(dims.nx), dims.K, fma), Ggt(vector<int>(dims.nu + dims.nx + 1), vector<int>(dims.ng), dims.K, fma)
        {
        };
        fatrop_memory_el<int> nu;
        fatrop_memory_el<int> nx;
        fatrop_memory_el<int> ng;
        fatrop_memory_matrix_bf_vec RSQrqt;
        fatrop_memory_matrix_bf_vec BAbt;
        fatrop_memory_matrix_bf_vec Ggt;
    };
} // namespace fatrop

#endif // FATROPOCPKKTINDLUCED