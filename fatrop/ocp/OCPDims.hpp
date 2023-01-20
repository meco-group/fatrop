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
#include "aux/FatropVector.hpp"
#include <vector>
using namespace std;
namespace fatrop
{
    /** \brief  this class contains the problem dimensions of a standard ocp*/
    struct OCPDims
    {
    public:
        OCPDims(const int K, const FatropVector<int> &nu, const FatropVector<int> &nx, const FatropVector<int> &ng, const FatropVector<int> &ng_ineq) : K(K), nu(nu), nx(nx), ng(ng), ng_ineq(ng_ineq),
                                                                                                                                                        n_u_tot(sum(nu)),
                                                                                                                                                        n_x_tot(sum(nx)),
                                                                                                                                                        n_b_tot(n_x_tot - nx.at(0)),
                                                                                                                                                        n_g_tot(sum(ng)),
                                                                                                                                                        n_g_ineq_tot(sum(ng_ineq)){};
        /// horizon length
        const int K;
        /// input vector size
        const FatropVector<int> nu;
        /// state vector size
        const FatropVector<int> nx;
        // number of stagewise equality constraints
        const FatropVector<int> ng;
        // number of stagewise inequality constraints
        const FatropVector<int> ng_ineq;
        const int n_u_tot;
        const int n_x_tot;
        const int n_b_tot;
        const int n_g_tot;
        const int n_g_ineq_tot;
    };
} // namespace fatrop
#endif // FATROP_OCPDIMS_INCLUDED