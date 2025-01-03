#ifndef __fatrop_unittest_random_matrix_hpp__
#define __fatrop_unittest_random_matrix_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/linear_algebra.hpp"
#include <random>
namespace fatrop::test
{
    // Function to generate a random double
    Scalar random(Scalar lower_bound = 0.0, Scalar upper_bound = 1.0)
    {
        // Static ensures the generator and distribution are initialized once
        static std::mt19937 gen(0); // Random number engine
        std::uniform_real_distribution<Scalar> dist(lower_bound, upper_bound);
        return dist(gen);
    }
    MatRealAllocated random_matrix(Index rows, Index cols, Scalar lower_bound = 0.0,
                                   Scalar upper_bound = 1.0)
    {
        MatRealAllocated matrix(rows, cols);
        for (Index i = 0; i < rows; ++i)
        {
            for (Index j = 0; j < cols; ++j)
            {
                matrix(i, j) = random(lower_bound, upper_bound);
            }
        }
        return matrix;
    }
    MatRealAllocated random_spd_matrix(Index m, Scalar lower_bound = 0.0, Scalar upper_bound = 1.0)
    {
        MatRealAllocated matrix = random_matrix(m, m, lower_bound, upper_bound);
        MatRealAllocated ret(m, m);
        syrk_ln(m, m, 1.0, matrix, 0, 0, matrix, 0, 0, 0.0, ret, 0, 0, ret, 0, 0);
        trtr_l(m, ret, 0, 0, ret, 0, 0);
        return ret;
    }
}
#endif  // __fatrop_unittest_random_matrix_hpp__