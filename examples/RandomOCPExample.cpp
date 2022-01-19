#include "DEBUG/RandomOCP.hpp"
#include "OCP/OCPLSRiccati.hpp"
#include "OCP/OCPAlg.hpp"
#include "AUX/FatropVector.hpp"
#include <vector>
#include <iostream>
#include "BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
using namespace std;
using namespace fatrop;
int main()
{
    int K = 10;
    int nu = 5;
    int nx = 10;
    int ng = 0;
    FatropVector<int> nu_ = vector<int>(K, nu);
    FatropVector<int> nx_ = vector<int>(K, nx);
    FatropVector<int> ng_ = vector<int>(K, ng);
    // ng_.at(K - 1) = nx;
    // ng_.at(0) = nx;

    MemoryAllocator fma;
    RefCountPtr<BFOCP> ocptemplatebasic =
        new RandomOCP(nu_, nx_, ng_, K);
    RefCountPtr<OCP> ocptempladapter = new BFOCPAdapter(ocptemplatebasic, fma);
    RefCountPtr<OCPLinearSolver> ocplsriccati = new OCPLSRiccati(ocptempladapter->GetOCPDims(), fma);
    FatropOCP ocpalg(ocptempladapter, ocplsriccati, fma);
    int N_opti_vars = sum(nu_ + nx_);
    int N_lags = sum(nx_) - nx_.at(0) + sum(ng_);
    FatropMemoryVecBF ux(N_opti_vars, 1, fma);
    FatropMemoryVecBF lags(N_lags, 1, fma);
    FatropMemoryVecBF lags2(N_lags, 1, fma);
    FatropMemoryVecBF ux2(N_opti_vars, 1, fma);
    fma.allocate();
    ocpalg.EvalHess(1.0, ux[0], ux[0], lags[0], lags[0]);
    ocpalg.EvalJac(ux[0], ux[0], lags[0]);
    ocplsriccati->computeSD(&ocpalg.ocpkktmemory_, 0.0, ux[0], lags[0]);
    ocpalg.ocpkktmemory_.BAbt[0].print();
    cout << Eig(ux[0]) << endl;
}