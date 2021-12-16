#ifndef SparseSolverInterfaceIncluded
#define SparseSolverInterfaceIncluded
#include "FatropSparse.hpp"
#include <vector>
#include "../AUX/FatropVector.hpp"
using namespace std;
namespace fatrop
{
    class SparseSolverInterface
    {
    public:
    // triplets are from lower part of matrix
        SparseSolverInterface(const int nnz, const int dim, const vector<triplet> &tripl) : nnz_(nnz), dim_(dim), ai(nnz), aj(nnz)
        {
            for (int i =0; i < nnz; i++)
            {
                ai.at(i) = tripl.at(i).ai ;
                aj.at(i) = tripl.at(i).aj ;
            }
        };
        virtual void preprocess() = 0;
        virtual void solve(const vector<triplet>& A, vector<double>& rhs) = 0;
    protected:
        const int nnz_;
        const int dim_;
        FatropVector<int> ai;
        FatropVector<int> aj;
    };
} //namespace fatrop

#endif // SparseSolverInterfaceIncluded