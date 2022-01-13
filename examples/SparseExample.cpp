#include "Fatrop.hpp"
#include "SPARSE/SparseOCP.hpp"
#include "SPARSE/InterfaceMUMPS.hpp"
#include "DEBUG/FatropDebugTools.hpp"
#include "AUX/FatropVector.hpp"
// #include <gperftools/profiler.h>
using namespace fatrop;
int main()
{
    // ProfilerStart("output_inside.prof"); //Start profiling section and save to file

    /// sparse ocp
    OCPDims dims;
    // next parameters give problem
    // dims.K = 3;
    // int nu = 10;
    // int nx = 9;
    // int ng = 4;
    // dims.nx.at(2) = 3*nx;
    dims.K = 100;
    int nu = 5;
    int nx = 10;
    int ng = 0;
    dims.nx = vector<int>(dims.K, nx);
    dims.nu = vector<int>(dims.K, nu);
    dims.ng = vector<int>(dims.K, ng);
    // dims.nu.at(dims.K - 1) = 0;
    dims.ng.at(dims.K - 1) = nx;
    dims.ng.at(0) = nx;
    // dims.ng.at(dims.K-1) = 0;
    // memory allocation
    MemoryAllocator fma;
    OCPKKT KKTocp(dims, fma);
    OCP_KKT_solver OCP_solver(dims, fma);

    int N_opti_vars = sum(dims.nu + dims.nx);
    int N_lags = sum(dims.nx) - dims.nx.at(0) + sum(dims.ng);
    FatropMemoryVecBF ux(N_opti_vars, 1, fma);
    FatropMemoryVecBF lags(N_lags, 1, fma);
    FatropMemoryVecBF lags2(N_lags, 1, fma);
    FatropMemoryVecBF ux2(N_opti_vars, 1, fma);
    fma.allocate();
    random_OCP(KKTocp, dims, 3);
    // KKTocp.BAbt[0].print();
    Sparse_OCP KOCP(dims, KKTocp);
    blasfeo_timer timer;
    // KOCP.KKT.print("matrix");
    double el;
    KOCP.fact_solve(ux[0], lags[0], el);
    cout << "solving using MUMPS" << endl;
    blasfeo_tic(&timer);
    KOCP.fact_solve(ux[0], lags[0], el);
    // lags[0].print();
    cout << "solving using fatrop" << endl;
    int N = 1000;
    blasfeo_tic(&timer);
    for (int i = 0; i < N; i++)
    {
        OCP_solver.fact_solve(&KKTocp, ux2[0], lags2[0]);
    }
    double el2 = blasfeo_toc(&timer) / N;
    // lags2[0].print();
    cout << "el time mumps " << el << endl;
    cout << "el time recursion " << el2 << endl;
    cout << "speedup " << el / el2 << endl;
    cout << "inf-norm difference MUMPS - Fatrop  primal " << (Eig(ux[0]) - Eig(ux2[0])).lpNorm<Eigen::Infinity>() << endl;
    cout << "inf-norm difference MUMPS - Fatrop  dual " << (Eig(lags[0]) - Eig(lags2[0])).lpNorm<Eigen::Infinity>() << endl;
    cout << "inf-norm   primal " << Eig(ux[0]).lpNorm<Eigen::Infinity>() << endl;
    cout << "inf-norm   dual " << Eig(lags[0]).lpNorm<Eigen::Infinity>() << endl;
    // ProfilerStop();
    return 0;
}