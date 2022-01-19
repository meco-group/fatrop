#include "Fatrop.hpp"
#include "DEBUG/FatropDebugTools.hpp"
#include "SPARSE/FatropSparse.hpp"
#include "vector"
#include "string"
#include <iostream>
#include "AUX/FatropVector.hpp"
#include <memory>
#include <OCP/BFOCPBasic.hpp>
#include <OCP/BFOCPAdapter.hpp>
using namespace std;
using namespace fatrop;

int main()
{
    // RefCountPtr<OCPTemplate> ocptempl = new OCPTemplateBasic(OCPTemplateBasic::from_shared_lib("f.so", )));
    // OCPDims dims;
    // // dims.K = 3;
    // // int nu = 2;
    // // int nx = 5;
    // // int ng = 1;
    // // dims.nx = vector<int>(dims.K, nx);
    // // dims.nu = vector<int>(dims.K, nu);
    // // dims.ng = vector<int>(dims.K, ng);
    // //     dims.K = 3;

    // dims.K = 5;
    // int nu = 10;
    // int nx = 9;
    // int ng = 5;
    // dims.nx = vector<int>(dims.K, nx);
    // dims.nu = vector<int>(dims.K, nu);
    // dims.ng = vector<int>(dims.K, ng);
    // dims.ng.at(2) = 0;
    // dims.ng.at(0) = 0;
    // // dims.ng.at(1) = 0;
    // dims.nx.at(2) = 25;
    // dims.nx.at(2) = 25;
    // // dims.ng.at(2) = 0;
    // // dims.ng.at(0) = 0;
    // // dims.nu.at(2) = 2*nu;
    // // dims.nu.at(20) = 0.5*nu;
    // dims.ng.at(dims.K-1) = 0;
    // dims.nu.at(dims.K-1) = 0;
    // // memory allocation
    // MemoryAllocator fma;
    // OCPKKTMemory KKT(dims, fma);
    // OCPKKTSolver OCP_solver(dims, fma);
    // int N_opti_vars = sum(dims.nu + dims.nx);
    // int N_lags = (dims.K - 1) * nx + sum(dims.ng);
    // FatropMemoryVecBF ux(N_opti_vars, 1, fma);
    // FatropMemoryVecBF lags(N_lags, 1, fma);
    // fma.allocate();
    // random_OCP(KKT, dims, 0);
    // blasfeo_timer timer;
    // blasfeo_tic(&timer);
    // for (int i = 0; i < 1; i++)
    // {
    //     OCP_solver.fact_solve(&KKT, ux[0], lags[0]);
    // }
    // double el = blasfeo_toc(&timer);
    // ux[0].print();
    // cout << "el time " << el << endl;
    return 0;
}