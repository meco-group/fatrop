/*
 * Fatrop - A fast trajectory optimization solver
 * Copyright (C) 2022, 2023 Lander Vanroye <lander.vanroye@kuleuven.be>
 *
 * This file is part of Fatrop.
 *
 * Fatrop is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fatrop is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Fatrop.  If not, see <http://www.gnu.org/licenses/>. */
#ifndef FATROP_OCPDIMS_INCLUDED
#define FATROP_OCPDIMS_INCLUDED
#include "auxiliary/FatropVector.hpp"
#include "templates/NLPAlg.hpp"
#include <vector>
namespace fatrop
{
    /** \brief  this class contains the problem dimensions of a standard ocp*/
    struct OCPDims : public NLPDims
    {
    public:
        OCPDims(const int K, const FatropVector<int> &nu, const FatropVector<int> &nx, const FatropVector<int> &ng, const FatropVector<int> &ng_ineq, const FatropVector<int> &n_stage_params, const int n_global_parameters) : NLPDims({sum(nx + nu), sum(ng + ng_ineq + nx) - nx.get(0), sum(ng_ineq)}),
                                                                                                                                                        K(K), nu(nu), nx(nx), ng(ng), ng_ineq(ng_ineq), n_stage_params(n_stage_params), n_global_params(n_global_parameters),
                                                                                                                                                        n_u_tot(sum(nu)),
                                                                                                                                                        n_x_tot(sum(nx)),
                                                                                                                                                        n_b_tot(n_x_tot - nx.get(0)),
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
        // number of stage parameters
        const FatropVector<int> n_stage_params;
        // number of global parameters
        const int n_global_params;
        const int n_u_tot;
        const int n_x_tot;
        const int n_b_tot;
        const int n_g_tot;
        const int n_g_ineq_tot;
    };
} // namespace fatrop
#endif // FATROP_OCPDIMS_INCLUDED