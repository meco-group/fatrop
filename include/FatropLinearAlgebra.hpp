/**
 *  @file FatropLinearAlgebra.hpp
 *  This file contains some interface classes for use in linear algebra. 
 *
*/
#ifndef FATROP_LA_INCLUDED
#define FATROP_LA_INCLUDED
namespace fatrop
{
    /** \brief interface class for matrix representations */
    class fatrop_matrix
    {
    public:
        /** \brief copy of matrix element */
        virtual double get_el(const int ai, const int aj) const = 0;
        /** \brief number of rows */
        virtual int nrows() const = 0;
        /** \brief number of cols */
        virtual int ncols() const = 0;
        void print()
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
    };
} // namespace fatrop

#endif // FATROP_LA_INCLUDED