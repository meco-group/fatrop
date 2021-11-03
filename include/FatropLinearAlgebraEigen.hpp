/**
 *  @file FatropLinearAlgebraEigen.hpp
 *  This file provides an interface for fatrop_matrix/vector to Eigen matrix resp. vecor, this is only used for debugging and testing purposes. 
 *  The code in this file is not optimized for efficiency!
*/
#ifndef FATROPLINEARALGEBRAEIGENINCLUDED
#define FATROPLINEARALGEBRAEIGENINCLUDED
#include <eigen3/Eigen/Dense>
#include "Fatrop.hpp"
namespace fatrop
{
    /** \brief convert fatrop_matrix to Eigen matrix, only used for testing/debugging*/
    Eigen::MatrixXd Eig(const fatrop_matrix &fm)
    {
        const int nrows = fm.nrows();
        const int ncols = fm.ncols();
        Eigen::MatrixXd res(nrows, ncols);
        for (int ai = 0; ai < nrows; ai++)
        {
            for (int aj = 0; aj < ncols; aj++)
            {
                res(ai, aj) = fm.get_el(ai, aj);
            }
        }
        return res;
    };
}; //namespace fatrop
#endif