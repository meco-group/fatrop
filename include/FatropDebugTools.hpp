#ifndef FATROPDEBUGTOOLSINCLUDED
#define FATROPDEBUGTOOLSINCLUDED
#include "FatropLinearAlgebraBlasfeo.hpp"
#include <random>
#include <vector>
void fill_matrix(const int m, const int n, MAT *A, const int ai, const int aj, const int seed = 0)
{
    std::default_random_engine e(seed);
    std::uniform_real_distribution<> dis(0, 1); // rage 0 - 1
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
        {
            {
                MATEL(A, ai + i, aj + j) = dis(e);
            }
        }
};
void fill_matrix(MAT *A, int seed = 0)
{
    fill_matrix(A->m, A->n, A, 0, 0, seed);
}
#endif