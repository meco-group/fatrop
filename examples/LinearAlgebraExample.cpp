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
    fatrop_memory_matrix_bf lower(10, 10, 1, fma);
    fatrop_memory_matrix_bf At(10, 10, 1, fma);
    fatrop_memory_matrix_bf A1(10, 10, 1, fma);
    fatrop_memory_matrix_bf L(10, 10, 1, fma);
    fatrop_memory_matrix_bf U(10, 10, 1, fma);
    fatrop_memory_matrix_bf test(10, 10, 1, fma);
    fatrop_memory_permutation_matrix Pl(10, 1, fma);
    fatrop_memory_permutation_matrix Pr(10, 1, fma);
    fma.allocate();
    fill_matrix((MAT *)A);
    At = Eig(Eig(A).transpose());
    GECP(10, 10, (MAT *)A, 0, 0, (MAT *)A1, 0, 0);
    cout << "A" << endl;
    A.print();
    cout << "At" << endl;
    At.print();
    int rank = 0;
    LU_FACT(10, 10, 10, rank, (MAT *)A, Pl.perm_vector(0), Pr.perm_vector(0));
    cout << "LU factorization " << endl;
    A.print();
    cout << "LUt factorization " << endl;
    LU_FACT_transposed(10, 10, 10, rank, (MAT *)At, Pl.perm_vector(0), Pr.perm_vector(0));
    At.print();
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            if (i == j)
            {
                MATEL((MAT *)L, i, i) = 1.0;
                MATEL((MAT *)U, i, i) = MATEL((MAT *)A, i, i);
            }
            else if (i > j)
            {
                // lower
                MATEL((MAT *)L, i, j) = MATEL((MAT *)A, i, j);
            }
            else if (i < j)
            {
                // upper
                MATEL((MAT *)U, i, j) = MATEL((MAT *)A, i, j);
            }
        }
    }
    cout << "L" << endl;
    L.print();
    cout << "U" << endl;
    U.print();
    cout << "A - Pl^T @ L @ U @ Pr" << endl;
    cout << Eig(A1) - Eig(Pl).transpose() * Eig(L) * Eig(U) * Eig(Pr) << endl;
    cout << "A - At^T" << Eig(Eig(A) - Eig(Eig(At).transpose())) <<endl;
    blasfeo_timer timer;
    blasfeo_tic(&timer);
    for(int i = 0; i < 1000; i++)
    {
        TRSM_RLNN(10, 10, 1.0, (MAT *)L, 0, 0, (MAT *)L, 0, 0, (MAT *)test, 0, 0);
    }
    double el = blasfeo_toc(&timer);
    test.print();
    cout << "el time " << el << endl;
    TRSM_RLNN(10, 10, 1.0, (MAT *)L, 0, 0, (MAT *)A, 0, 0, (MAT *)test, 0, 0);
    cout << Eig(Eig(test) - Eig(A)*Eig(Eig(L).inverse())) << endl;

    return 0;
}