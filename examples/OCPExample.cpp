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
    fma.allocate();
    random_OCP(KKT, dims, 0);
    blasfeo_dvec *dummy;
    KKT.RSQrqt[0].print();
    OCP_solver.fact_solve(&KKT, dummy, dummy);
    return 0;
}