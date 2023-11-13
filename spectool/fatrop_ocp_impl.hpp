#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "casadi_utilities.hpp"
#include "casadi_fe.hpp"
#include "ustage.hpp"
#include "ocp.hpp"
#include "ustage_quanitities.hpp"
#include "ocp/OCPAbstract.hpp"
#include "fatrop_ustage_eval_casadi.hpp"
namespace fatrop
{
    namespace spectool
    {
        namespace cs = casadi;
        class FatropOcpImpl : public OCPAbstract
        {
        public:
            FatropOcpImpl(const Ocp &ocp, const cs::Dict &opts)
            {
                horizon_length_ = 0;
                n_global_parameters_ = cs::MX::veccat(ocp.get_global_parameters()).size1();
                CasadiJitCache eval_cache;
                for (const auto &ustage : ocp.get_ustages())
                {
                    ustages_.push_back(std::make_shared<FatropuStageEvalCasadi>(uStageQuantities::create(ustage.get_internal()), opts, eval_cache));
                    for (int i = 1; i < ustage.K(); i++)
                        ustages_.push_back(ustages_.back());
                    horizon_length_ += ustage.K();
                }
            }
            fatrop_int get_nxk(const fatrop_int k) const override { return ustages_[k]->get_nx(); };
            fatrop_int get_nuk(const fatrop_int k) const override { return ustages_[k]->get_nu(); };
            fatrop_int get_ngk(const fatrop_int k) const override { return ustages_[k]->get_ng(); };
            fatrop_int get_n_stage_params_k(const fatrop_int k) const override { return ustages_[k]->get_n_stage_params(); };
            fatrop_int get_n_global_params() const override { return n_global_parameters_; };
            fatrop_int get_default_stage_paramsk(double *stage_params, const fatrop_int k) const override { return 0; };
            fatrop_int get_default_global_params(double *global_params) const override { return 0; };
            fatrop_int get_ng_ineq_k(const fatrop_int k) const override { return ustages_[k]->get_ng_ineq(); };
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
                return ustages_[k]->eval_BAbt(states_kp1, inputs_k, states_k, stage_params_k, global_params, res);
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
                return ustages_[k]->eval_RSQrqt(objective_scale, inputs_k, states_k, lam_dyn_k, lam_eq_k, lam_eq_ineq_k, stage_params_k, global_params, res);
            };
            fatrop_int eval_Ggtk(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res,
                const fatrop_int k) override
            {
                return ustages_[k]->eval_Ggt(inputs_k, states_k, stage_params_k, global_params, res);
            };
            fatrop_int eval_Ggt_ineqk(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res,
                const fatrop_int k) override
            {
                return ustages_[k]->eval_Ggt_ineq(inputs_k, states_k, stage_params_k, global_params, res);
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
                return ustages_[k]->eval_b(states_kp1, inputs_k, states_k, stage_params_k, global_params, res);
            };
            fatrop_int eval_gk(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res,
                const fatrop_int k) override
            {
                return ustages_[k]->eval_g(inputs_k, states_k, stage_params_k, global_params, res);
            };
            fatrop_int eval_gineqk(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res,
                const fatrop_int k) override
            {
                return ustages_[k]->eval_gineq(inputs_k, states_k, stage_params_k, global_params, res);
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
                return ustages_[k]->eval_rq(objective_scale, inputs_k, states_k, stage_params_k, global_params, res);
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
                return ustages_[k]->eval_L(objective_scale, inputs_k, states_k, stage_params_k, global_params, res);
            };
            fatrop_int get_boundsk(double *lower, double *upper, const fatrop_int k) const override
            {
                return ustages_[k]->get_bounds(lower, upper);
            };
            fatrop_int get_initial_xk(double *xk, const fatrop_int k) const override
            {
                // initialization is done higher up
                for (int i = 0; i < get_nxk(k); i++)
                {
                    xk[i] = 0;
                }
                return 0;
            };
            fatrop_int get_initial_uk(double *uk, const fatrop_int k) const override
            {
                // initialization is done higher up
                for (int i = 0; i < get_nuk(k); i++)
                {
                    uk[i] = 0;
                }
                return 0;
            };

        private:
            std::vector<std::shared_ptr<FatropuStageEvalAbstract>> ustages_;
            int horizon_length_ = 0;
            int n_global_parameters_ = 0;
        };
        // implementation of OCPAbstract, given an OCP
    } // namespace spectool
} // namespace fatrop