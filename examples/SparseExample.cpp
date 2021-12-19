#include "Fatrop.hpp"
#include "SPARSE/SparseOCP.hpp"
#include "SPARSE/InterfaceMUMPS.hpp"
#include "DEBUG/FatropDebugTools.hpp"
#include "AUX/FatropVector.hpp"
using namespace fatrop;
int main()
{
    // int N = 2;
    // KKT_matrix KKT;
    // var_sp x = KKT.get_variable(N);
    // KKT.set_hess_block(Eig(Id(N, N)), x, x);
    // x->set_grad(vector<double>(N, 0.0));
    // KKT.set_equation(Eig(Id(1, N)) * x, vector<double>(1, 1.0));
    // vector<triplet> tripvec;
    // KKT.get_triplets(tripvec);
    // KKT.print("matrix");
    // InterfaceMUMPS interf(tripvec.size(), KKT.get_size(), tripvec);
    // interf.preprocess();
    // vector<double> rhs(KKT.get_rhs());
    // interf.solve(tripvec, rhs);

    /// sparse ocp
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
    OCP_KKT KKTocp(dims, fma);
    int N_opti_vars = sum(dims.nu + dims.nx);
    int N_lags = (dims.K - 1) * nx + sum(dims.ng);
    fatrop_memory_vector_bf ux(N_opti_vars, 1, fma);
    fatrop_memory_vector_bf lags(N_lags, 1, fma);
    fma.allocate();
    random_OCP(KKTocp, dims, 0);
    Sparse_OCP KOCP(dims, KKTocp);
    KOCP.KKT.print("matrix");
    KOCP.fact_solve(ux[0], lags[0]);
    ux[0].print();

    vector<triplet> ocptripl;
    KOCP.KKT.get_triplets(ocptripl);
    InterfaceMUMPS interfo(ocptripl.size(), KOCP.KKT.get_size(), ocptripl);
    interfo.preprocess();
    vector<double> rhso(KOCP.KKT.get_rhs());
    interfo.solve(ocptripl, rhso);

    return 0;
}