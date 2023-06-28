#include "auxiliary/LinearAlgebra.hpp"
namespace fatrop
{
    void FatropMat::print()
    {
        int n_rows = nrows();
        int n_cols = ncols();
        for (int ai = 0; ai < n_rows; ai++)
        {
            for (int aj = 0; aj < n_cols; aj++)
            {
                printf("%9.5f ", get_el(ai, aj));
            }
            printf("\n");
        }
    }
    void FatropVec::print()
    {
        int n_el = nels();
        for (int ai = 0; ai < n_el; ai++)
        {
            printf("%9.5f\n", get_el(ai));
        }
    }
}
