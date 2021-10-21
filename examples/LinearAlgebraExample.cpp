#include <iostream>
#include <memory>
#include "Fatrop.hpp"
using namespace fatrop;
using namespace std;
class test_container: public fatrop_memory_allocator{
};
int main(){
    fatrop_memory_allocator fma;
    fatrop_memory_el<int> test(5, fma);
    fma.allocate();
    std::cout << *((int*) test) << std::endl;
    return 0;
}