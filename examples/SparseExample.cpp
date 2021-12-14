#include "Fatrop.hpp"
#include "FatropSparse.hpp"
#include "SparseSolvers/InterfaceMUMPS.hpp"
#include "FatropDebugTools.hpp"
using namespace fatrop;
int main()
{
    int N = 2;
    KKT_matrix KKT;
    var_sp x = KKT.get_variable(N);
    KKT.set_hess_block(Eig(Id(N, N)), x, x);
    x->set_grad(vector<double>(N, 0.0));
    KKT.set_equation(Eig(Id(1, N)) * x, vector<double>(1, 1.0));
    vector<triplet> tripvec;
    KKT.get_triplets(tripvec);
    KKT.print("matrix");
    InterfaceMUMPS interf(tripvec.size(), KKT.get_size(), tripvec);
    interf.preprocess();
    vector<double> rhs(KKT.get_rhs());
    interf.solve(tripvec, rhs);

    /// sparse ocp
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
    OCP_KKT KKTocp(dims, fma);
    fma.allocate();
    random_OCP(KKTocp, dims, 0);
    KKT_matrix KOCP(Sparse_OCP(dims, KKTocp));
    vector<triplet> ocptripl;
    KOCP.get_triplets(ocptripl);
    InterfaceMUMPS interfo(ocptripl.size(),KOCP.get_size(), ocptripl);
    interfo.preprocess();
    vector<double> rhso(KOCP.get_rhs());
    interfo.solve(ocptripl, rhso);
    return 0;
}