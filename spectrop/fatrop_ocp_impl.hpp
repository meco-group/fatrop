#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "casadi_utilities.hpp"
#include "casadi_fe.hpp"
#include "ustage.hpp"
#include "ocp.hpp"
#include "ustage_quanitities.hpp"
#include "ocp/OCPAbstract.hpp"
namespace fatrop
{
    namespace spectrop
    {
        namespace cs = casadi;
        struct FatropuStageEval
        {
            CasadiFEWrap RSQrqt;
            CasadiFEWrap BAbt;
            CasadiFEWrap L;
            CasadiFEWrap b; // actually provided as dynamics
            CasadiFEWrap g_equality;
            CasadiFEWrap g_inequality;
            CasadiFEWrap rq;
            CasadiFEWrap Ggt_equality;
            CasadiFEWrap Ggt_inequality;
            std::vector<double> Lb;
            std::vector<double> Ub;
            int K;
            int nu;
            int nx;
            int np_stage;
            int np_global;
            int ng_eq;
            int ng_ineq;
            static FatropuStageEval create(const uStage &ustage, const cs::Dict& opts);
        };
        class FatropOcpImpl : public OCPAbstract
        {
        public:
            FatropOcpImpl(const Ocp &ocp, const cs::Dict &opts)
            {
                horizon_length_ = 0;
                for (const auto &ustage : ocp.get_ustages())
                {
                    ustages_.push_back(std::make_shared<FatropuStageEval>(FatropuStageEval::create(ustage, opts)));
                    for (int i = 1; i < ustage.K(); i++)
                        ustages_.push_back(ustages_.back());
                    horizon_length_ += ustage.K();
                }
            }
            fatrop_int get_nxk(const fatrop_int k) const override { return ustages_[k]->nx; };
            fatrop_int get_nuk(const fatrop_int k) const override { return ustages_[k]->nu; };
            fatrop_int get_ngk(const fatrop_int k) const override { return ustages_[k]->ng_eq; };
            fatrop_int get_n_stage_params_k(const fatrop_int k) const override { return ustages_[k]->np_stage; };
            fatrop_int get_n_global_params() const override { return ustages_[0]->np_global; };
            fatrop_int get_default_stage_paramsk(double *stage_params, const fatrop_int k) const override { return 0; };
            fatrop_int get_default_global_params(double *global_params) const override { return 0; };
            fatrop_int get_ng_ineq_k(const fatrop_int k) const override { return ustages_[k]->ng_ineq; };
            fatrop_int get_horizon_length() const override { return horizon_length_; };
            fatrop_int eval_BAbtk(
                const double *states_kp1,
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res,
                const fatrop_int k) override
            {

                const double *arg[5];
                arg[0] = states_k;
                arg[1] = states_kp1;
                arg[2] = inputs_k;
                arg[3] = stage_params_k;
                arg[4] = global_params;
                ustages_[k]->BAbt.eval(arg, res);
                return 0;
            };
            fatrop_int eval_RSQrqtk(
                const double *objective_scale,
                const double *inputs_k,
                const double *states_k,
                const double *lam_dyn_k,
                const double *lam_eq_k,
                const double *lam_eq_ineq_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res,
                const fatrop_int k) override
            {
                const double *arg[7];
                arg[0] = inputs_k;
                arg[1] = states_k;
                arg[2] = lam_dyn_k;
                arg[3] = lam_eq_k;
                arg[4] = lam_eq_ineq_k;
                arg[5] = stage_params_k;
                arg[6] = global_params;
                ustages_[k]->RSQrqt.eval(arg, res);
                return 0;
            };
            fatrop_int eval_Ggtk(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res,
                const fatrop_int k) override
            {
                const double *arg[4];
                arg[0] = inputs_k;
                arg[1] = states_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                ustages_[k]->Ggt_equality.eval(arg, res);
                return 0;
            };
            fatrop_int eval_Ggt_ineqk(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res,
                const fatrop_int k) override
            {
                const double *arg[4];
                arg[0] = inputs_k;
                arg[1] = states_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                ustages_[k]->Ggt_inequality.eval(arg, res);
                return 0;
            };
            fatrop_int eval_bk(
                const double *states_kp1,
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res,
                const fatrop_int k) override
            {
                const double *arg[5];
                arg[0] = states_k;
                arg[1] = states_kp1;
                arg[2] = inputs_k;
                arg[3] = stage_params_k;
                arg[4] = global_params;
                ustages_[k]->b.eval(arg, res);
                return 0;
            };
            fatrop_int eval_gk(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res,
                const fatrop_int k) override
            {
                const double *arg[4];
                arg[0] = inputs_k;
                arg[1] = states_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                ustages_[k]->g_equality.eval(arg, res);

                return 0;
            };
            fatrop_int eval_gineqk(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res,
                const fatrop_int k) override
            {
                const double *arg[4];
                arg[0] = inputs_k;
                arg[1] = states_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                ustages_[k]->g_inequality.eval(arg, res);

                return 0;
            };
            fatrop_int eval_rqk(
                const double *objective_scale,
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res,
                const fatrop_int k) override
            {
                const double *arg[4];
                arg[0] = inputs_k;
                arg[1] = states_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                ustages_[k]->rq.eval(arg, res);
                return 0;
            };
            fatrop_int eval_Lk(
                const double *objective_scale,
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res,
                const fatrop_int k) override
            {
                const double *arg[4];
                arg[0] = inputs_k;
                arg[1] = states_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                ustages_[k]->L.eval(arg, res);
                return 0;
            };
            fatrop_int get_boundsk(double *lower, double *upper, const fatrop_int k) const override
            {
                for (int i = 0; i < ustages_[k]->ng_ineq; i++)
                {
                    lower[i] = ustages_[k]->Lb[i];
                    upper[i] = ustages_[k]->Ub[i];
                }
                return 0;
            };
            fatrop_int get_initial_xk(double *xk, const fatrop_int k) const override
            {
                // initialization is done higher up
                for (int i = 0; i < ustages_[k]->nx; i++)
                {
                    xk[i] = 0;
                }
                return 0;
            };
            fatrop_int get_initial_uk(double *uk, const fatrop_int k) const override
            {
                // initialization is done higher up
                for (int i = 0; i < ustages_[k]->nu; i++)
                {
                    uk[i] = 0;
                }
                return 0;
            };

        private:
            std::vector<std::shared_ptr<FatropuStageEval>> ustages_;
            int horizon_length_ = 0;
        };
        // implementation of OCPAbstract, given an OCP
    } // namespace spectrop
} // namespace fatrop