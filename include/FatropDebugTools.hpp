#ifndef FATROPDEBUGTOOLSINCLUDED
#define FATROPDEBUGTOOLSINCLUDED
#include "FatropLinearAlgebraBlasfeo.hpp"
#include <random>
#include <vector>
#include "FatropLinearAlgebraEigen.hpp"
using namespace fatrop;
Eig random_matrix(const int m, const int n, const int seed = 0)
{
    Eig res(m, n);
    std::default_random_engine e(seed);
    std::uniform_real_distribution<> dis(0, 1); // rage 0 - 1
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
        {
            {
                res(i, j) = dis(e);
            }
        }
    return res;
}
OCP_KKT random_OCP(OCP_KKT& KKT, OCP_dims &dims, int seed = 0)
{
    int K = dims.K;
    for (int k = 0; k < K - 1; k++)
    {
        int nu = dims.nu.at(k);
        int nx = dims.nx.at(k);
        int ng = dims.ng.at(k);
        Eig rand1 = random_matrix(nu + nx, nu + nx, seed + 0);
        KKT.RSQrqt[k].block(0, 0, nu + nx, nu + nx) = Eig(rand1.transpose() * rand1);
        KKT.RSQrqt[k].block(nu + nx, 0, 1, nu + nx) = random_matrix(1, nu + nx, seed + 1);
        KKT.BAbt[k] = random_matrix(nu + nx + 1, nx, seed + 2);
        if (ng > 0)
        {
            KKT.Ggt[k] = random_matrix(nu + nx + 1, ng, seed + 3);
        }
    }
    // k = K-1
    int nu = dims.nu.at(K - 1);
    int nx = dims.nx.at(K - 1);
    int ng = dims.ng.at(K - 1);
    Eig rand1 = random_matrix(nx, nx, seed + 3);
    KKT.RSQrqt[K - 1].block(nu, nu, nx, nx) = Eig(rand1.transpose() * rand1);
    KKT.RSQrqt[K-1].block(nu + nx, nu, 1, nx) = random_matrix(1, nx, seed + 3);
    KKT.Ggt[K-1].block(nu, 0, nx+1, ng) = random_matrix(nx+1, ng, seed + 4);
    return KKT;
}
#endif