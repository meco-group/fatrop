#include "Fatrop.hpp"
#include "FatropDebugTools.hpp"
#include "FatropSparse.hpp"
#include "vector"
#include "string"
#include <iostream>
#include "FatropVector.hpp"
using namespace std;
using namespace fatrop;

int main()
{
    OCP_dims dims;
    dims.K = 10;
    int nu = 5;
    int nx = 9;
    int ng = 2;
    dims.nx = vector<int>(dims.K, nx);
    dims.nu = vector<int>(dims.K, nu);
    dims.ng = vector<int>(dims.K, ng);
    // memory allocation
    fatrop_memory_allocator fma;
    OCP_KKT KKT(dims, fma);
    OCP_KKT_solver OCP_solver(dims, fma);
    int N_opti_vars = dims.K*nx+(dims.K-1)*nu;
    int N_lags = (dims.K-1)*nx + sum(dims.ng);
    fatrop_memory_vector_bf ux(N_opti_vars, 1, fma);
    fatrop_memory_vector_bf lags(N_lags, 1, fma);
    fma.allocate();
    random_OCP(KKT, dims, 0);
    blasfeo_dvec *dummy;
    KKT.RSQrqt[0].print();
    OCP_solver.fact_solve(&KKT, (VEC*) ux, (VEC*) lags);
    return 0;
}