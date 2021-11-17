#include "Fatrop.hpp"
#include "FatropDebugTools.hpp"
#include "FatropSparse.hpp"
#include "vector"
using namespace std;
using namespace fatrop;
int main(){
    // problem dimensions
    OCP_dims dims;
    dims.K = 3;
    dims.nx = vector<int>(dims.K, 3); 
    dims.nu = vector<int>(dims.K, 3); 
    dims.ng = vector<int>(dims.K, 1); 
    // memory allocation
    fatrop_memory_allocator fma;
    OCP_KKT KKT(dims, fma);
    fma.allocate();
    random_OCP(KKT, dims, 0);
    Sparse_OCP(dims, KKT).print("matrix");
    return 0;
}