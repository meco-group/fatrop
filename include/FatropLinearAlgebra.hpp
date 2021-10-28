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
        /** \brief acces to matrix element */
        virtual double& at(const int ai, const int aj) const = 0;
        /** \brief number of rows */
        virtual int nrows() const = 0;
        /** \brief number of cols */
        virtual int ncols() const = 0;
    };
} // namespace fatrop

#endif // FATROP_LA_INCLUDED