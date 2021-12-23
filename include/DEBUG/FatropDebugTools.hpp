#ifndef FATROPDEBUGTOOLSINCLUDED
#define FATROPDEBUGTOOLSINCLUDED
#include "../BLASFEO_WRAPPER/FatropLinearAlgebraBlasfeo.hpp"
#include <random>
#include <vector>
#include "FatropLinearAlgebraEigen.hpp"
using namespace fatrop;
Eig random_matrix(const int m, const int n, const int seed = 0)
{
    Eig res(m, n);
    std::default_random_engine e(seed);
    std::uniform_real_distribution<> dis(0, 1); // rage 0 - 1
    for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
        {
            {
                res(i, j) = dis(e);
            }
        }
    return res;
}
OCP_KKT random_OCP(OCP_KKT &KKT, OCP_dims &dims, int seed = 0)
{
    int K = dims.K;
    for (int k = 0; k < K - 1; k++)
    {
        int nu = dims.nu.at(k);
        int nx = dims.nx.at(k);
        int nxp1 = dims.nx.at(k+1);
        int ng = dims.ng.at(k);
        KKT.RSQrqt[k] = random_matrix(nu + nx + 1, nu + nx, seed + 0);
        Eig rand1 = Eig(KKT.RSQrqt[k].block(0, 0, nu + nx, nu + nx));
        KKT.RSQrqt[k].block(0, 0, nu + nx, nu + nx) = Eig(rand1.transpose() * rand1);
        KKT.BAbt[k] = random_matrix(nu + nx + 1, nxp1, seed + 1);
        if (ng > 0)
        {
            KKT.Ggt[k].block(nu, 0, nx + 1, ng) = random_matrix(nx + 1, ng, seed + 2);
        }
    }
    // k = K-1
    int nu = dims.nu.at(K - 1);
    int nx = dims.nx.at(K - 1);
    int ng = dims.ng.at(K - 1);
    KKT.RSQrqt[K - 1].block(nu, nu, nx + 1, nx) = random_matrix(nx + 1, nx, seed + 0);
    Eig rand1 = Eig(KKT.RSQrqt[K - 1].block(nu, nu, nx, nx));
    KKT.RSQrqt[K - 1].block(nu, nu, nx, nx) = Eig(rand1.transpose() * rand1);
    KKT.Ggt[K - 1].block(nu, 0, nx + 1, ng) = random_matrix(nx + 1, ng, seed + 2);
    return KKT;
}
#endif