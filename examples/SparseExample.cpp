#include "Fatrop.hpp"
#include "FatropSparse.hpp"
#include "SparseSolvers/InterfaceMUMPS.hpp"
using namespace fatrop;
int main()
{
    int N = 2;
    KKT_matrix KKT;
    var_sp x = KKT.get_variable(N);
    KKT.set_hess_block(Eig(Id(N, N)), x, x);
    x ->set_grad(vector<double>(N,0.0));
    KKT.set_equation(Eig(Id(1,N))*x, vector<double>(1,1.0));
    vector<triplet> tripvec;
    KKT.get_triplets(tripvec);
    KKT.print("matrix");
    InterfaceMUMPS interf(tripvec.size(), KKT.get_size(), KKT);
    interf.preprocess();
    return 0;
}