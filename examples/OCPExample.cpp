#include "Fatrop.hpp"
#include "DEBUG/FatropDebugTools.hpp"
#include "SPARSE/FatropSparse.hpp"
#include "vector"
#include "string"
#include <iostream>
#include "AUX/FatropVector.hpp"
using namespace std;
using namespace fatrop;

int main()
{
    OCP_dims dims;
    dims.K = 3;
    int nu = 2;
    int nx = 2;
    int ng = 1;
    dims.nx = vector<int>(dims.K, nx);
    dims.nu = vector<int>(dims.K, nu);
    dims.ng = vector<int>(dims.K, ng);
    // memory allocation
    fatrop_memory_allocator fma;
    OCP_KKT KKT(dims, fma);
    OCP_KKT_solver OCP_solver(dims, fma);
    int N_opti_vars = sum(dims.nu+dims.nx);
    int N_lags = (dims.K-1)*nx + sum(dims.ng);
    fatrop_memory_vector_bf ux(N_opti_vars, 1, fma);
    fatrop_memory_vector_bf lags(N_lags, 1, fma);
    fma.allocate();
    random_OCP(KKT, dims, 0);
    cout << "RSQrqt \n" << endl;
    KKT.RSQrqt[0].print();
    cout << "BAbt \n" << endl;
    KKT.BAbt[0].print();
    cout << "Ggt \n" << endl;
    KKT.Ggt[0].print();
    cout << "---------------- \n" << endl;
    OCP_solver.fact_solve(&KKT, ux[0], lags[0]);
    ux[0].print();
    return 0;
}