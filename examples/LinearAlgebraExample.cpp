#include <iostream>
#include <memory>
#include "Fatrop.hpp"
#include "FatropLinearAlgebraEigen.hpp"
#include "FatropDebugTools.hpp"
using namespace fatrop;
using namespace std;
int main()
{
    fatrop_memory_allocator fma;
    fatrop_memory_matrix_bf A(10, 10, 1, fma);
    fatrop_memory_matrix_bf lower(vector<int>(1, 10), vector<int>(1, 10), 1, fma);
    fatrop_memory_matrix_bf At(10, 10, 1, fma);
    fatrop_memory_matrix_bf A1(10, 10, 1, fma);
    fatrop_memory_matrix_bf L(10, 10, 1, fma);
    fatrop_memory_matrix_bf U(10, 10, 1, fma);
    fatrop_memory_matrix_bf test(10, 10, 1, fma);
    fatrop_memory_permutation_matrix Pl(10, 1, fma);
    fatrop_memory_permutation_matrix Pr(10, 1, fma);
    fatrop_memory_el<int> testfatropmemel(5, vector<int>(5, 420), fma);
    fma.allocate();
    A[0] = random_matrix(10, 10);
    At[0] = Eig(Eig(A[0]).transpose());
    GECP(10, 10, (MAT *)A[0], 0, 0, (MAT *)A1[0], 0, 0);
    cout << "A" << endl;
    A[0].print();
    cout << "At" << endl;
    At[0].print();
    int rank = 0;
    LU_FACT(10, 10, 10, rank, (MAT *)A[0], ((PMAT *)Pl), ((PMAT *)Pr));
    cout << "LU factorization " << endl;
    A[0].print();
    cout << "LUt factorization " << endl;
    LU_FACT_transposed(10, 10, 10, rank, (MAT *)At[0], ((PMAT *)Pl), ((PMAT *)Pr));
    At[0].print();
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            if (i == j)
            {
                MATEL((MAT *)L[0], i, i) = 1.0;
                MATEL((MAT *)U[0], i, i) = MATEL((MAT *)A[0], i, i);
            }
            else if (i > j)
            {
                // lower
                MATEL((MAT *)L[0], i, j) = MATEL((MAT *)A[0], i, j);
            }
            else if (i < j)
            {
                // upper
                MATEL((MAT *)U[0], i, j) = MATEL((MAT *)A[0], i, j);
            }
        }
    }
    cout << "L" << endl;
    L[0].print();
    cout << "U" << endl;
    U[0].print();
    cout << "A - Pl^T @ L @ U @ Pr" << endl;
    cout << Eig(A1[0]) - Eig(Pl).transpose() * Eig(L[0]) * Eig(U[0]) * Eig(Pr) << endl;
    cout << "A - At^T" << Eig(Eig(A[0]) - Eig(Eig(At[0]).transpose())) << endl;
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    for (int i = 0; i < 1000; i++)
    {
        TRSM_RLNN(10, 10, 1.0, (MAT *)L[0], 0, 0, (MAT *)L[0], 0, 0, (MAT *)test[0], 0, 0);
    }
    double el = blasfeo_toc(&timer);
    test[0].print();
    cout << "el time " << el << endl;
    TRSM_RLNN(10, 10, 1.0, (MAT *)L[0], 0, 0, (MAT *)A[0], 0, 0, (MAT *)test[0], 0, 0);
    cout << Eig(Eig(test[0]) - Eig(A[0]) * Eig(Eig(L[0]).inverse())) << endl;
    cout << ((int *)testfatropmemel)[0] << endl;

    cout << Eig(test[0]) << std::endl;

    return 0;
}