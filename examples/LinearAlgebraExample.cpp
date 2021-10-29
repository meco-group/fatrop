#include <iostream>
#include <memory>
#include "Fatrop.hpp"
using namespace fatrop;
using namespace std;
class test_container : public fatrop_memory_allocator
{
};
int main()
{
    fatrop_memory_allocator fma;
    fatrop_memory_el<int> test(5, fma);
    fatrop_memory_matrix test2(5, 5, 1, fma);
    fatrop_memory_permutation_matrix P(4, 1, fma);
    fma.allocate();
    int *P_data = P.perm_vector(0);
    P_data[0] = 2;
    P_data[1] = 2;
    P_data[2] = 3;
    P_data[3] = 3;
    P.print();
    std::cout << *((int *)test) << std::endl;
    return 0;
}