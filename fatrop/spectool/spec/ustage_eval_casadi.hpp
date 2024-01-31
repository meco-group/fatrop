#pragma once
extern "C"
{
#include <blasfeo.h>
}
#include "fatrop/ocp/UStageEvalAbstract.hpp"
#include "fatrop/spectool/function_evaluation/casadi_fe.hpp"
#include "fatrop/spectool/spec/ustage_quantities.hpp"
namespace fatrop
{
    namespace spectool
    {
        class FatropuStageEvalCasadi : public UStageEvalAbstract
        {
        public:
            FatropuStageEvalCasadi(const uStageQuantities &sq, const cs::Dict &opts, CasadiJitCache &eval_cache);
            virtual int nu(const uStageQuantities &sq);
            virtual int nx(const uStageQuantities &sq);
            virtual int np_stage(const uStageQuantities &sq);
            virtual int np_global(const uStageQuantities &sq);
            virtual int ng_eq(const uStageQuantities &sq);
            virtual int ng_ineq(const uStageQuantities &sq);
            virtual int nxp1(const uStageQuantities &sq);

        public:
            virtual int get_nx() const;
            virtual int get_nu() const;
            virtual int get_ng() const;
            virtual int get_n_stage_params() const;
            virtual int get_ng_ineq() const;
            virtual int eval_BAbt(
                const double *states_kp1,
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res);
            virtual int eval_RSQrqt(
                const double *objective_scale,
                const double *inputs_k,
                const double *states_k,
                const double *lam_dyn_k,
                const double *lam_eq_k,
                const double *lam_eq_ineq_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res);
            virtual int eval_Ggt(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res);
            virtual int eval_Ggt_ineq(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res);
            virtual int eval_b(
                const double *states_kp1,
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res);
            virtual int eval_g(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res);
            virtual int eval_gineq(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res);
            virtual int eval_rq(
                const double *objective_scale,
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res);
            virtual int eval_L(
                const double *objective_scale,
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res);
            virtual int get_bounds(double *lower, double *upper) const;
            CasadiFEWrap RSQrqt_;
            CasadiFEWrap BAbt_;
            CasadiFEWrap L_;
            CasadiFEWrap b_; // actually provided as dynamics
            CasadiFEWrap g_equality_;
            CasadiFEWrap g_inequality_;
            CasadiFEWrap rq_;
            CasadiFEWrap Ggt_equality_;
            CasadiFEWrap Ggt_inequality_;
            std::vector<double> Lb_;
            std::vector<double> Ub_;
            int K_;
            int nu_;
            int nx_;
            int nxp1_;
            int np_stage_;
            int np_global_;
            int ng_eq_;
            int ng_ineq_;
        };
    }
}