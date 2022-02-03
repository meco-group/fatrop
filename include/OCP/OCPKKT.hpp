/**
 * @file OCPKKT.hpp
 * @author your name (you@domain.com)
 * @brief this file contains classes to represent a KKT matrix by blasfeo submatrices
 * @version 0.1
 * @date 2021-11-10
 *
 * @copyright Copyright (c) 2021
 *
 */
#ifndef FATROP_OCP_KKT_INCLUDED
#define FATROP_OCP_KKT_INCLUDED
#include "../BLASFEO_WRAPPER/LinearAlgebraBlasfeo.hpp"
#include "OCPDims.hpp"
#include "../AUX/Aux.hpp"
#include "../AUX/FatropVector.hpp"
using namespace std;
namespace fatrop
{
    /** \brief this class contains all information to represent the KKT system of an equality constrained OCP*/
    class OCPKKTMemory
    {
    public:
        OCPKKTMemory(const OCPDims &dims) : K(dims.K),
                                            nu(dims.nu),
                                            nx(dims.nx),
                                            ng(dims.ng),
                                            RSQrqt(dims.nu + dims.nx + 1, dims.nu + dims.nx, dims.K),
                                            BAbt(dims.nu + dims.nx + 1, rotate(dims.nx, 1), dims.K),
                                            Ggt(dims.nu + dims.nx + 1, dims.ng, dims.K),
                                            Ggt_ineq(dims.nu + dims.nx + 1, dims.ng_ineq, dims.K),
                                            aux(dims){};

        OCPKKTMemory(const OCPKKTMemory &cpy) = delete;
        int K;
        FatropVector<int> nu;
        FatropVector<int> nx;
        FatropVector<int> ng;
        /// small-scale Hessian
        FatropMemoryMatBF RSQrqt;
        /// small-scale Jacobian dynamics
        FatropMemoryMatBF BAbt;
        /// small-scale Jacobian stagewise eq constraints
        FatropMemoryMatBF Ggt;
        /// small-scale Jacobian stagewise ineq constraints
        FatropMemoryMatBF Ggt_ineq;
        class OCPAux
        {
        public:
            OCPAux(const OCPDims &dims) : ux_offs(offsets(dims.nx + dims.nu)),
                                          g_offs(offsets(dims.ng)),
                                          dyn_eq_offs(offsets(rotate(dims.nx, 1)) + sum(dims.ng)),
                                          g_ineq_offs(offsets(rotate(dims.ng_ineq, 1)) + (sum(dims.nx) - dims.nx.at(0) + sum(dims.ng))),
                                          max_nu(max(dims.nu)), max_nx(max(dims.nx)),
                                          max_ng(max(dims.ng)),
                                          max_ngineq(max(dims.ng_ineq)){};
            /// offset arrays are used for efficiency
            const FatropVector<int> ux_offs;
            /// offset arrays are used for efficiency
            const FatropVector<int> g_offs;
            const FatropVector<int> dyn_eq_offs;
            const FatropVector<int> g_ineq_offs;
            int max_nu;
            int max_nx;
            int max_ng;
            int max_ngineq;
        };
        OCPAux aux;
    };

} // namespace fatrop
#endif // FATROP_OCP_KKT_INCLUDED