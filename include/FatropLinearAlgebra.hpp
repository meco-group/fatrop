/**
 *  @file FatropLinearAlgebra.hpp
 *  This file contains some interface classes for use in linear algebra. 
 *
*/
#ifndef FATROP_LA_INCLUDED
#define FATROP_LA_INCLUDED
namespace fatrop
{
    struct matrix_ind
    {
        int ai;
        int aj;
    };
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
    // special matrices
    class eye : public fatrop_matrix
    {
    public:
        /** \brief constructor */
        eye(int dim) : dim_(dim){};
        /** \brief copy of matrix element */
        double get_el(const int ai, const int aj) const
        {
            if (ai == aj)
            {
                return 1.0;
            }
            else
            {
                return 0.0;
            }
        };
        /** \brief number of rows */
        int nrows() const { return dim_; };
        /** \brief number of cols */
        int ncols() const { return dim_; };

    private:
        int dim_;
    };
    /** \brief interface class for matrix representations */
    class fatrop_vector
    {
    public:
        /** \brief copy of matrix element */
        virtual double get_el(const int ai) const = 0;
        /** \brief number of elements */
        virtual int nels() const = 0;
        void print()
        {
            int n_el = nels();
            for (int ai = 0; ai < n_el; ai++)
            {
                printf("%9.5f ", get_el(ai));
                printf("\n");
            }
        }
    };
} // namespace fatrop

#endif // FATROP_LA_INCLUDED