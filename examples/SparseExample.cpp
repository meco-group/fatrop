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
    int nx = 4;
    int ng = 2;
    dims.nx = vector<int>(dims.K, nx);
    dims.nu = vector<int>(dims.K, nu);
    dims.ng = vector<int>(dims.K, 0);
    dims.ng.at(0) = ng;
    dims.nu.at(dims.K-1) = 0;
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
    KKTocp.BAbt[0].print();
    Sparse_OCP KOCP(dims, KKTocp);
    blasfeo_timer timer;
    cout << "solving using MUMPS" << endl;
    // KOCP.KKT.print("matrix");
    blasfeo_tic(&timer);
    KOCP.fact_solve(ux[0], lags[0]);
    lags[0].print();
    double el = blasfeo_toc(&timer);
    cout << "solving using fatrop" << endl;
    blasfeo_tic(&timer);
    OCP_solver.fact_solve(&KKTocp, ux2[0], lags[0]);
    lags[0].print();
    double el2 = blasfeo_toc(&timer);
    cout << "el time mumps " << el << endl;
    cout << "el time recursion " << el2 << endl;
    cout << "inf-norm difference MUMPS - Fatrop " << (Eig(ux[0]) - Eig(ux2[0])).lpNorm<Eigen::Infinity>() << endl;
    return 0;
}