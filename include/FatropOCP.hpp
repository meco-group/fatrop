/**
 * @file FatropOCP.hpp
 * @author your name (you@domain.com)
 * @brief this file contains necessary the OCP specific code 
 * @version 0.1
 * @date 2021-11-10
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef FATROP_OCP_INCLUDED
#define FATROP_OCP_INCLUDED
#include "FatropMemory.hpp"
#include "FatropVector.hpp"
#include <vector>
using namespace std;
namespace fatrop
{
    /** \brief  this class contains the problem dimensions of a standard ocp*/
    struct OCP_dims
    {
    public:
        /// horizon length
        int K;
        /// input vector size
        FatropVector<int> nu;
        /// state vector size
        FatropVector<int> nx;
        // number of stagewise equality constraints
        FatropVector<int> ng;
    };
} // namespace fatrop
#endif //FATROPOCPINCLUDED