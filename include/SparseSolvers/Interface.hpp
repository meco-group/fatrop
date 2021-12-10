#ifndef SparseSolverInterfaceIncluded
#define SparseSolverInterfaceIncluded
#include "../FatropSparse.hpp"
#include <vector>
using namespace std;
namespace fatrop
{
    class SparseSolverInterface
    {
    public:
        SparseSolverInterface(const int nnz, const int dim, KKT_matrix &KKT_m) : nnz_(nnz), dim_(dim), ai(dim), aj(dim)
        {
            vector<triplet> tripl;
            KKT_m.get_triplets(tripl);
            for (int i; i < dim; i++)
            {
                ai.at(i) = tripl.at(i).ai;
                aj.at(i) = tripl.at(i).aj;
            }
        };
        virtual void preprocess() = 0;
        virtual void solve(const vector<triplet> A, const vector<double> rhs) = 0;
    private:
        const int nnz_;
        const int dim_;
        vector<int> ai;
        vector<int> aj;
    };
} //namespace fatrop

#endif // SparseSolverInterfaceIncluded