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
using namespace std;
namespace fatrop
{
    class OCP_KKT
    {
    public:
        OCP_KKT(const fatrop_OCP_dims &dims, fatrop_memory_allocator &fma) : dims(dims), RSQrqt(dims.K, fma), BAbt(dims.K, fma), Ggt(dims.K, fma)
        {
            int K = dims.K;
            int nu = dims.nu;
            int nx = dims.nx;
            RSQrqt.set_dimensions(vector<int>(K, nu+nx+1), vector<int>(K, nu+nx));
            BAbt.set_dimensions(vector<int>(K, nu+nx+1), vector<int>(K, nu+nx));
            Ggt.set_dimensions(vector<int>(K, nu+nx+1), vector<int>(dims.ng));
        };

    private:
        fatrop_OCP_dims dims;
        fatrop_memory_matrix_bf_vec RSQrqt;
        fatrop_memory_matrix_bf_vec BAbt;
        fatrop_memory_matrix_bf_vec Ggt;
    };
} // namespace fatrop

#endif // FATROPOCPKKTINDLUCED