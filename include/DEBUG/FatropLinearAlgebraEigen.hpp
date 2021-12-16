/**
 *  @file FatropLinearAlgebraEigen.hpp
 *  This file provides an interface for fatrop_matrix/vector to Eigen matrix resp. vecor, this is only used for debugging and testing purposes. 
 *  The code in this file is not optimized for efficiency!
*/
#ifndef FATROP_LINEAR_ALGEBRA_EIGEN_INCLUDED
#define FATROP_LINEAR_ALGEBRA_EIGEN_INCLUDED
#include <eigen3/Eigen/Dense>
#include "../Fatrop.hpp"
namespace fatrop
{

    /** \brief convert fatrop_matrix to Eigen matrix, only used for testing/debugging*/
    class Eig : public Eigen::MatrixXd, public fatrop_matrix
    {
    public:
        Eig(const fatrop_matrix &fm) : Eigen::MatrixXd(fm.nrows(), fm.ncols())
        {
            const int nrows = fm.nrows();
            const int ncols = fm.ncols();
            for (int ai = 0; ai < nrows; ai++)
            {
                for (int aj = 0; aj < ncols; aj++)
                {
                    Eigen::MatrixXd::operator()(ai, aj) = fm.get_el(ai, aj);
                }
            }
        };
        Eig(const Eigen::MatrixXd &Eigenmat) : Eigen::MatrixXd(Eigenmat){};
        Eig(const int m, const int n) : Eigen::MatrixXd(m, n){};
        /** \brief a normal Eig reprensents a dense matrix */
        virtual bool iszero(const int ai, const int aj) const { return false; };
        /** \brief copy of matrix element */
        double get_el(const int ai, const int aj) const { return Eigen::MatrixXd::operator()(ai, aj); };
        /** \brief number of rows */
        int nrows() const { return Eigen::MatrixXd::rows(); };
        /** \brief number of cols */
        int ncols() const { return Eigen::MatrixXd::cols(); };
    };
    // /** \brief Eig matrix with diagonal sparse structure */
    // class EigD : public Eig
    // {
    //     public:
    //     EigD(const fatrop_matrix &fm) : Eig(fm){};
    //     bool iszero(const int ai, const int aj) const { return (ai == aj) ? false : true; cout << "here" << endl;}
    // };
}; //namespace fatrop
#endif