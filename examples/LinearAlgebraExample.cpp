#include <iostream>
#include <memory>
#include "Fatrop.hpp"
#include "debug/LinearAlgebraEigen.hpp"
#include "debug/FatropDebugTools.hpp"
#include <limits>
using namespace fatrop;
using namespace std;
int main()
{
    blasfeo_dmat sA;
    blasfeo_dmat sB;
    blasfeo_dmat sC;
    blasfeo_dmat sD;
    blasfeo_allocate_dmat(12, 12, &sA);
    blasfeo_allocate_dmat(12, 12, &sB);
    blasfeo_allocate_dmat(12, 12, &sC);
    blasfeo_allocate_dmat(12, 12, &sD);
    GESE(12, 12, 0.0, &sA, 0, 0);
    GESE(12, 12, 0.0, &sB, 0, 0);
    GESE(12, 12, 0.0, &sC, 0, 0); // 12x12
    GESE(12, 12, 0.0, &sD, 0, 0);
    // D <- A@B^T + C
    GEMM_NT(11, 11, 11, 1.0, &sA, 0, 0, &sB, 1, 0, 1.0, &sC, 0, 0, &sD, 0, 0);
    blasfeo_free_dmat(&sA);
    blasfeo_free_dmat(&sB);
    blasfeo_free_dmat(&sC);
    blasfeo_free_dmat(&sD);

    // FatropMemoryMatBF A(100, 100, 1);
    // FatropMemoryMatBF lower(vector<int>(1, 100), vector<int>(1, 100), 1);
    // FatropMemoryMatBF At(100, 100, 1);
    // FatropMemoryMatBF A1(100, 100, 1);
    // FatropMemoryMatBF L(100, 100, 1);
    // FatropMemoryMatBF U(100, 100, 1);
    // FatropMemoryMatBF test(100, 100, 1);
    // MemoryPermMat Pl(100, 1);
    // MemoryPermMat Pr(100, 1);
    // vector<int> testfatropmemel(5, 420);
    // A[0] = random_matrix(100, 100);
    // At[0] = Eig(Eig(A[0]).transpose());
    // GECP(100, 100, (MAT *)A[0], 0, 0, (MAT *)A1[0], 0, 0);
    // // cout << "A" << endl;
    // // A[0].print();
    // // cout << "At" << endl;
    // // At[0].print();
    // int rank = 0;
    // LU_FACT(100, 100, 100, rank, (MAT *)A[0], ((PMAT *)Pl), ((PMAT *)Pr));
    // // cout << "LU factorization " << endl;
    // // A[0].print();
    // // cout << "LUt factorization " << endl;
    // LU_FACT_transposed(100, 100, 100, rank, (MAT *)At[0], ((PMAT *)Pl), ((PMAT *)Pr));
    // // At[0].print();
    // for (int i = 0; i < 100; i++)
    // {
    //     for (int j = 0; j < 100; j++)
    //     {
    //         if (i == j)
    //         {
    //             MATEL((MAT *)L[0], i, i) = 1.0;
    //             MATEL((MAT *)U[0], i, i) = MATEL((MAT *)A[0], i, i);
    //         }
    //         else if (i > j)
    //         {
    //             // lower
    //             MATEL((MAT *)L[0], i, j) = MATEL((MAT *)A[0], i, j);
    //         }
    //         else if (i < j)
    //         {
    //             // upper
    //             MATEL((MAT *)U[0], i, j) = MATEL((MAT *)A[0], i, j);
    //         }
    //     }
    // }
    // // cout << "L" << endl;
    // // L[0].print();
    // // cout << "U" << endl;
    // // U[0].print();
    // // cout << "A - Pl^T @ L @ U @ Pr" << endl;
    // // cout << Eig(A1[0]) - Eig(Pl).transpose() * Eig(L[0]) * Eig(U[0]) * Eig(Pr) << endl;
    // // cout << "A - At^T" << Eig(Eig(A[0]) - Eig(Eig(At[0]).transpose())) << endl;
    // blasfeo_timer timer;
    // MAT* L0 = (MAT *)L[0];
    // MAT* test0 = (MAT *)test[0];
    // blasfeo_tic(&timer);
    // for (int i = 0; i < 1000; i++)
    // {
    //     TRSM_RLNN(100, 100, 1.0,L0 , 0, 0,L0, 0, 0, test0, 0, 0);
    // }
    // double el = blasfeo_toc(&timer);
    // // test[0].print();
    // cout << "el time " << el << " s" << endl;
    // TRSM_RLNN(100, 100, 1.0, (MAT *)L[0], 0, 0, (MAT *)A[0], 0, 0, (MAT *)test[0], 0, 0);
    // // cout << Eig(Eig(test[0]) - Eig(A[0]) * Eig(Eig(L[0]).inverse())) << endl;
    // // cout << ((int *)testfatropmemel)[0] << endl;

    // // cout << Eig(test[0]) << std::endl;

    return 0;
}