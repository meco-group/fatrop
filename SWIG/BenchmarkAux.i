%module BenchmarkAux
    %{
    #define SWIG_FILE_WITH_INIT
    #include "DEBUG/RandomOCP.hpp"
    #include "AUX/SmartPtr.hpp"
    #include <vector>
    RefCountPtr<RandomOCP> GenerateRandom(int nu, int nx, int ng, int K)
    {
    FatropVector<int> nu_ = vector<int>(K, nu);
    FatropVector<int> nx_ = vector<int>(K, nx);
    FatropVector<int> ng_ = vector<int>(K, ng);
        return new RandomOCP(nu_, nx_, ng_, K);
    }
    %}

    RefCountPtr<RandomOCP> GenerateRandom(int nu, int nx, int ng, int K);
