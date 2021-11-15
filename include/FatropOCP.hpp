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
#ifndef FATROPOCPINCLUDED
#define FATROPOCPINCLUDED
#include "FatropMemory.hpp"
#include <vector>
using namespace std;
namespace fatrop
{
    /** \brief  this class contains the problem dimensions of a standard ocp*/
    struct fatrop_OCP_dims
    {
    public:
        /// horizon length
        int K;
        /// input vector size
        vector<int> nu;
        /// state vector size
        vector<int> nx;
        // number of stagewise equality constraints
        vector<int> ng;
    };
} // namespace fatrop
#endif //FATROPOCPINCLUDED