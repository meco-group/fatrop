#pragma once
#include "fatrop-casadi-problem.hpp"
#include <ocp/OCPAbstract.hpp>
namespace fatrop
{
    namespace specification
    {
        class FatropCasadiSolver : public OCPAbstract
        {
        public:
            fatrop_int get_nxk(const fatrop_int k) const override
            {
                return evaluation_quantities_ptr_[k]->nx;
            };
            fatrop_int get_nuk(const fatrop_int k) const override
            {
                return evaluation_quantities_ptr_[k]->nu;
            };
            fatrop_int get_ngk(const fatrop_int k) const override
            {
                return evaluation_quantities_ptr_[k]->ng_eq;
            };
            fatrop_int get_ng_ineq_k(const fatrop_int k) const override
            {
                return evaluation_quantities_ptr_[k]->ng_ineq;
            };
            fatrop_int get_n_global_params() const
            {
                return np_global_;
            };
            fatrop_int get_n_stage_params_k(const fatrop_int k) const override
            {
                return evaluation_quantities_ptr_[k]->np_stage;
            };
            fatrop_int get_horizon_length() const override
            {
                return horizon_length_;
            };
            // implementation of the OCPAbstract methods
            int eval_BAbtk(const double *states_kp1,
                           const double *inputs_k,
                           const double *states_k,
                           const double *stage_params_k,
                           const double *global_params,
                           MAT *res,
                           const int k) override
            {
                const double *arg[5];
                arg[0] = states_k;
                arg[1] = states_kp1;
                arg[2] = inputs_k;
                arg[3] = stage_params_k;
                arg[4] = global_params;
                evaluation_quantities_ptr_[k]->BAbt(arg, res);
                return 0;
            };
            int eval_RSQrqtk(const double *objective_scale,
                             const double *inputs_k,
                             const double *states_k,
                             const double *lam_dyn_k,
                             const double *lam_eq_k,
                             const double *lam_ineq_k,
                             const double *stage_params_k,
                             const double *global_params,
                             MAT *res,
                             const int k) override
            {
                const double *arg[7];
                arg[0] = states_k;
                arg[1] = inputs_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                arg[4] = lam_dyn_k;
                arg[5] = lam_eq_k;
                arg[6] = lam_ineq_k;
                evaluation_quantities_ptr_[k]->RSQrq(arg, res);
                return 0;
            }
            int eval_Ggtk(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res,
                const int k) override
            {
                const double *arg[4];
                arg[0] = states_k;
                arg[1] = inputs_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                evaluation_quantities_ptr_[k]->Ggt_equality(arg, res);
                return 0;
            }
            int eval_Ggt_ineqk(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res,
                const int k)
            {
                const double *arg[4];
                arg[0] = states_k;
                arg[1] = inputs_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                evaluation_quantities_ptr_[k]->Ggt_inequality(arg, res);
                return 0;
            }
            int eval_bk(
                const double *states_kp1,
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res,
                const int k) override
            {
                const double *arg[5];
                arg[0] = states_k;
                arg[1] = states_kp1;
                arg[2] = inputs_k;
                arg[3] = stage_params_k;
                arg[4] = global_params;
                evaluation_quantities_ptr_[k]->b(arg, res);
                return 0;
            }
            int eval_gk(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res,
                const int k) override
            {
                const double *arg[4];
                arg[0] = states_k;
                arg[1] = inputs_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                evaluation_quantities_ptr_[k]->g_equality(arg, res);
                return 0;
            }
            int eval_gineqk(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res,
                const int k) override
            {
                const double *arg[4];
                arg[0] = states_k;
                arg[1] = inputs_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                evaluation_quantities_ptr_[k]->g_inequality(arg, res);
                return 0;
            }
            int eval_rqk(
                const double *objective_scale,
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res,
                const int k) override
            {
                const double *arg[4];
                arg[0] = states_k;
                arg[1] = inputs_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                evaluation_quantities_ptr_[k]->rq(arg, res);
                return 0;
            }

            int eval_Lk(
                const double *objective_scale,
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res,
                const int k) override
            {
                const double *arg[4];
                arg[0] = states_k;
                arg[1] = inputs_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                evaluation_quantities_ptr_[k]->L(arg, res);
                return 0;
            }
            int get_initial_xk(double *xk, const int k) const override
            {
                int nx = this->get_nxk(k);
                for (int i = 0; i < nx; i++)
                    xk[i] = x_initial[k][i];
                return 0;
            };
            int get_initial_uk(double *uk, const int k) const override
            {
                int nu = this->get_nuk(k);
                for (int i = 0; i < nu; i++)
                    uk[i] = u_initial[k][i];
                return 0;
            };
            int set_initial_xk(double *xk, const int k)
            {
                int nx = this->get_nxk(k);
                for (int i = 0; i < nx; i++)
                    x_initial[k][i] = xk[i];
                return 0;
            };
            int set_initial_uk(double *uk, const int k)
            {
                int nu = this->get_nuk(k);
                for (int i = 0; i < nu; i++)
                    u_initial[k][i] = uk[i];
                return 0;
            };
            int get_boundsk(double *lower, double *upper, const int k) const override
            {
                int ngineq = this->get_ng_ineq_k(k);
                const double *lower_ineq = evaluation_quantities_ptr_[k]->lb.data();
                const double *upper_ineq = evaluation_quantities_ptr_[k]->ub.data();

                for (int j = 0; j < ngineq; j++)
                {
                    lower[j] = lower_ineq[j];
                    upper[j] = upper_ineq[j];
                }
                return 0;
            };
            int get_default_stage_paramsk(double *stage_params_res, const int k) const override
            {
                int np = this->get_n_stage_params_k(k);
                const double *p_stage_default = evaluation_quantities_ptr_[k]->p_stage_default.data();
                for (int i = 0; i < np; i++)
                    stage_params_res[i] = p_stage_default[i];
                return 0;
            }
            int get_default_global_params(double *global_params_res) const override
            {
                int np = this->get_n_global_params();
                const double *p_global_default = p_global_default_.data();
                for (int i = 0; i < np; i++)
                    global_params_res[i] = p_global_default[i];
                return 0;
            }

        private:
            typedef std::vector<double> vd;
            struct evaluation_quantities
            {
                int nx, nu, np_stage, ng_eq, ng_ineq;
                eval_bf L;
                eval_bf RSQrq;
                eval_bf rq;
                eval_bf g_equality;
                eval_bf g_inequality;
                eval_bf Ggt_equality;
                eval_bf Ggt_inequality;
                eval_bf BAbt;
                eval_bf b;
                vd lb;
                vd ub;
                vd p_stage_default;
            };
            std::vector<evaluation_quantities> evaluation_quantities_;
            std::vector<evaluation_quantities *> evaluation_quantities_ptr_;
            vd p_global_default_;
            int np_global_;
            int horizon_length_;
            std::vector<vd> x_initial;
            std::vector<vd> u_initial;
        };
    } // namespace specification
} // namespace fatrop