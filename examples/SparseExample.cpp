#include "Fatrop.hpp"
#include "SPARSE/SparseOCP.hpp"
#include "SPARSE/InterfaceMUMPS.hpp"
#include "DEBUG/FatropDebugTools.hpp"
#include "AUX/FatropVector.hpp"
using namespace fatrop;
int main()
{
    /// sparse ocp
    OCP_dims dims;
    dims.K = 3;
    int nu = 2;
    int nx = 5;
    int ng = 1;
    dims.nx = vector<int>(dims.K, nx);
    dims.nu = vector<int>(dims.K, nu);
    dims.ng = vector<int>(dims.K, ng);
    // memory allocation
    fatrop_memory_allocator fma;
    OCP_KKT KKTocp(dims, fma);
    OCP_KKT_solver OCP_solver(dims, fma);

    int N_opti_vars = sum(dims.nu + dims.nx);
    int N_lags = (dims.K - 1) * nx + sum(dims.ng);
    fatrop_memory_vector_bf ux(N_opti_vars, 1, fma);
    fatrop_memory_vector_bf lags(N_lags, 1, fma);
    fatrop_memory_vector_bf ux2(N_opti_vars, 1, fma);
    fma.allocate();
    random_OCP(KKTocp, dims, 0);
    Sparse_OCP KOCP(dims, KKTocp);
    blasfeo_timer timer;
    // for (int i = 0; i < 1000; i++)
    // {
    //     KOCP.fact_solve(ux[0], lags[0]);
    // }
    blasfeo_tic(&timer);
    for (int i = 0; i < 1; i++)
    {
        KOCP.fact_solve(ux[0], lags[0]);
    }
    double el = blasfeo_toc(&timer);
    ux[0].print();
    // for (int i = 0; i < 1000; i++)
    // {
    //     OCP_solver.fact_solve(&KKTocp, ux2[0], lags[0]);
    // }
    blasfeo_tic(&timer);
    for (int i = 0; i < 1; i++)
    {
        OCP_solver.fact_solve(&KKTocp, ux2[0], lags[0]);
    }
    double el2 = blasfeo_toc(&timer);
    cout << "-----\n";
    ux2[0].print();
    cout << "el time mumps " << el << endl;
    cout << "el time recursion " << el2 << endl;
    return 0;
}