#ifndef SparseSolverInterfaceIncluded
#define SparseSolverInterfaceIncluded
#include "../FatropSparse.hpp"
#include <vector>
#include "FatropVector.hpp"
using namespace std;
namespace fatrop
{
    class SparseSolverInterface
    {
    public:
    // triplets are from lower part of matrix
        SparseSolverInterface(const int nnz, const int dim, KKT_matrix &KKT_m) : nnz_(nnz), dim_(dim), ai(nnz), aj(nnz)
        {
            vector<triplet> tripl;
            KKT_m.get_triplets(tripl);
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