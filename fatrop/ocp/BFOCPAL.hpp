// takes a BFOCP and transforms it to another BFOCP that represents it augmented Lagrangian subproblem
// it's iself a BFOPC as well such that it can be solved using the regular methods

#ifndef BFOCPALINCLUDED
#define BFOCPALINCLUDED
#include "BFOCP.hpp"
#include <memory>
#include <algorithm>
#include <cmath>
namespace fatrop
{
    class BFOCPAL : public BFOCP
    {
    public:
        BFOCPAL(shared_ptr<BFOCP> &ocp, double penalty) : ocp_(ocp), K(ocp_->get_horizon_length()),
                                                          no_ineqs(TransformRange<int>(0, K, [&ocp](int k)
                                                                                       { return ocp->get_ng_ineq_k(k); })),
                                                          nu(TransformRange<int>(0, K, [&ocp](int k)
                                                                                 { return ocp->get_nxk(k); })),
                                                          nx(TransformRange<int>(0, K, [&ocp](int k)
                                                                                 { return ocp->get_nuk(k); })),
                                                          ineqs_offsets(offsets(no_ineqs)),
                                                          ineq_lagsL(sum(no_ineqs), 1),
                                                          ineq_lagsU(sum(no_ineqs), 1),
                                                          lower_bounds(sum(no_ineqs), 1),
                                                          upper_bounds(sum(no_ineqs), 1),
                                                          penalty(penalty),
                                                          tmpmat(max(nu + nx) +1, max(no_ineqs), 1),
                                                          tmpviolation(max(no_ineqs), 1),
                                                          gradvec(max(no_ineqs), 1),
                                                          lagsupdated(max(no_ineqs), 1){};
        int get_nxk(const int k) const;
        int get_nuk(const int k) const;
        int get_ngk(const int k) const;
        int get_n_stage_params_k(const int k) const;
        int get_n_global_parmas() const;
        int get_ng_ineq_k(const int k) const;
        int get_ng_ineq_k_AL(const int k) const;
        int get_horizon_length() const;
        int eval_BAbtk(
            const double *states_kp1,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params_k,
            MAT *res,
            const int k);
        int eval_RSQrqtk(
            const double *objective_scale,
            const double *inputs_k,
            const double *states_k,
            const double *lam_dyn_k,
            const double *lam_eq_k,
            const double *lam_eq_ineq_k,
            const double *stage_params_k,
            const double *global_params_k,
            MAT *res,
            const int k);
        int eval_Ggtk(
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params_k,
            MAT *res,
            const int k);
        int eval_Ggt_ineqk(
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params_k,
            MAT *res,
            const int k);
        int eval_Ggt_ineqk_AL(
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params_k,
            MAT *res,
            const int k);
        int eval_bk(
            const double *states_kp1,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params_k,
            double *res,
            const int k);
        int eval_gk(
            const double *states_k,
            const double *inputs_k,
            const double *stage_params_k,
            const double *global_params_k,
            double *res,
            const int k);
        int eval_gineqk(
            const double *states_k,
            const double *inputs_k,
            const double *stage_params_k,
            const double *global_params_k,
            double *res,
            const int k);
        int eval_gineqk_AL(
            const double *states_k,
            const double *inputs_k,
            const double *stage_params_k,
            const double *global_params_k,
            double *res,
            const int k);
        int eval_rqk(
            const double *objective_scale,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params_k,
            double *res,
            const int k);
        int eval_Lk(
            const double *objective_scale,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params_k,
            double *res,
            const int k);
        shared_ptr<BFOCP> ocp_;
        const int K;
        // vector with number of ineqs each stage
        FatropVector<int> no_ineqs;
        // vector with number of inputs
        FatropVector<int> nu;
        // vector with number of states
        FatropVector<int> nx;
        // vector with ineqs offsets
        FatropVector<int> ineqs_offsets;
        // vector with inequality lags
        FatropMemoryVecBF ineq_lagsL;
        // vector with inequality lags
        FatropMemoryVecBF ineq_lagsU;
        // vector with inequality lags
        FatropMemoryVecBF lower_bounds;
        // vector with inequality lags
        FatropMemoryVecBF upper_bounds;
        // penalty parameter
        double penalty = 0.0;
        // blasfeo matrix for temp results
        FatropMemoryMatBF tmpmat;
        // blasfeo matrix for temp results
        FatropMemoryVecBF tmpviolation;
        // blasfeo matrix for temp results
        FatropMemoryVecBF gradvec;
        // blasfeo matrix for temp results
        FatropMemoryVecBF lagsupdated;
    };
} // namespace fatrop
#endif // BFOCPALINCLUDED
