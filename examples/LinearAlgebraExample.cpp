#include <iostream>
#include <memory>
#include "Fatrop.hpp"
#include "FatropLinearAlgebraEigen.hpp"
using namespace fatrop;
using namespace std;
class test_container : public fatrop_memory_allocator
{
};
int main()
{
    fatrop_memory_allocator fma;
    fatrop_memory_el<int> test(5, fma);
    fatrop_memory_matrix_bf test2(4, 4, 1, fma);
    fatrop_memory_matrix_bf test69(4, 4, 1, fma);
    fatrop_memory_permutation_matrix P(4, 1, fma);
    fma.allocate();
    int *P_data = P.perm_vector(0);
    P_data[0] = 2;
    P_data[1] = 2;
    P_data[2] = 3;
    P_data[3] = 3;
    P.print();
    fatrop_matrix_bf test3(4,4,0,0);
    test2 = P;
    blasfeo_print_dmat(4,4,(blasfeo_dmat*) test2, 0,0);
    std::cout << *((int *)test) << std::endl;
    test69 = eye(4);
    P.PM(4, (blasfeo_dmat*) test69);
    blasfeo_print_dmat(4,4,(blasfeo_dmat*) test69, 0,0);
    cout << Eig(test69) << std::endl;
    cout << Eig(test69.block(0,0,2,2)) << std::endl;
    return 0;
}