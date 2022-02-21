/**
 *  @file LinearAlgebra.hpp
 *  This file contains some interface classes for use in linear algebra.
 *
 */
#ifndef FATROP_LA_INCLUDED
#define FATROP_LA_INCLUDED
#include <iostream>
namespace fatrop
{
    struct MatrixInd
    {
        int ai;
        int aj;
    };
    /** \brief Interface class for matrix representations */
    class FatropMat
    {
    public:
        /** \brief Copy of matrix element */
        virtual double get_el(const int ai, const int aj) const = 0;
        /** \brief Number of rows */
        virtual int nrows() const = 0;
        /** \brief Number of cols */
        virtual int ncols() const = 0;
        void print();
    };
    // special matrices
    /** \brief Identity matrix */
    class eye : public FatropMat
    {
    public:
        /** \brief Constructor */
        eye(int dim) : dim_(dim){};
        /** \brief Copy of matrix element */
        inline double get_el(const int ai, const int aj) const
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
        /** \brief Number of rows */
        int nrows() const { return dim_; };
        /** \brief Number of cols */
        int ncols() const { return dim_; };

    private:
        int dim_;
    };
    /** \brief Interface class for matrix representations */
    class FatropVec
    {
    public:
        /** \brief Copy of matrix element */
        virtual double get_el(const int ai) const = 0;
        /** \brief Number of elements */
        virtual int nels() const = 0;
        void print();
    };
} // namespace fatrop

#endif // FATROP_LA_INCLUDED