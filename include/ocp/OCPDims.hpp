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
#ifndef FATROP_OCPDIMS_INCLUDED
#define FATROP_OCPDIMS_INCLUDED
#include "../aux/FatropMemory.hpp"
#include "../aux/FatropVector.hpp"
#include <vector>
using namespace std;
namespace fatrop
{
    /** \brief  this class contains the problem dimensions of a standard ocp*/
    struct OCPDims
    {
    public:
        OCPDims(){};
        OCPDims(const int K, const FatropVector<int>&nu,const FatropVector<int>& nx, const FatropVector<int>& ng):K(K), nu(nu), nx(nx), ng(ng){};
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
#endif // FATROP_OCPDIMS_INCLUDED