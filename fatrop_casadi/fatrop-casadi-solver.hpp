#pragma once
#include "fatrop-casadi-problem.hpp"
#include "function-evaluation.hpp"
#include <ocp/OCPAbstract.hpp>
#include <unordered_map>
namespace fatrop
{
    namespace fatrop_casadi
    {
        class FatropCasadiSolver : public OCPAbstract
        {
        public:
            FatropCasadiSolver(const FatropCasadiProblem& prob)
            {
                horizon_length_ = prob.size();
                np_global_ = prob[0] -> dims.np_global;
                u_initial.reserve(horizon_length_);
                x_initial.reserve(horizon_length_ + 1);
                std::unordered_map<std::shared_ptr<MicroStageInternal>, std::shared_ptr<evaluation_quantities> > evaluation_quantities_map;
                for (int k = 0; k < horizon_length_; k++)
                {
                    // check if prob[k] is already in the map
                    if (evaluation_quantities_map.find(prob[k]) == evaluation_quantities_map.end())
                    {
                        // if not, create a new entry
                        evaluation_quantities_map[prob[k]] = std::make_shared<evaluation_quantities>(prob[k]);
                    }
                    evaluation_quantities_ptr_.push_back(evaluation_quantities_map[prob[k]]);
                    // zero initialize x and u 
                    int nu = evaluation_quantities_map[prob[k]]->nu;
                    int nx = evaluation_quantities_map[prob[k]]->nx;
                    u_initial.push_back(std::vector<double>(nu, 0.0));
                    x_initial.push_back(std::vector<double>(nx, 0.0));
                }
            };
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
                evaluation_quantities_ptr_[k]->BAbt->operator()(arg, res);
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
                evaluation_quantities_ptr_[k]->RSQrq->operator()(arg, res);
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
                evaluation_quantities_ptr_[k]->Ggt_equality->operator()(arg, res);
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
                evaluation_quantities_ptr_[k]->Ggt_inequality->operator()(arg, res);
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
                evaluation_quantities_ptr_[k]->b->operator()(arg, res);
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
                evaluation_quantities_ptr_[k]->g_equality->operator()(arg, res);
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
                evaluation_quantities_ptr_[k]->g_inequality->operator()(arg, res);
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
                evaluation_quantities_ptr_[k]->rq->operator()(arg, res);
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
                evaluation_quantities_ptr_[k]->L->operator()(arg, res);
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
                evaluation_quantities(const MicroStage& stagequantities)
                {
                    nx = stagequantities -> dims.nx;
                    int nxp1 = stagequantities -> dims.nxp1;
                    nu = stagequantities -> dims.nu;
                    np_stage = stagequantities -> dims.np_stage;
                    ng_eq = stagequantities -> dims.ng_equality;
                    ng_ineq = stagequantities -> dims.ng_inequality;
                    L = EvalBF::create(stagequantities -> L);
                    RSQrq = EvalBF::create(stagequantities -> RSQrq);
                    rq = EvalBF::create(stagequantities -> rq);
                    if(ng_eq>0)
                    g_equality = EvalBF::create(stagequantities -> g_equality);
                    if(ng_ineq>0)
                    g_inequality = EvalBF::create(stagequantities -> g_inequality);
                    if(ng_eq>0)
                    Ggt_equality = EvalBF::create(stagequantities -> Ggt_equality);
                    if(ng_ineq>0)
                    Ggt_inequality = EvalBF::create(stagequantities -> Ggt_inequality);
                    if(nxp1>0)
                    BAbt = EvalBF::create(stagequantities -> BAbt);
                    if(nxp1>0)
                    b = EvalBF::create(stagequantities -> b);
                    lb = stagequantities -> Lb;
                    ub = stagequantities -> Ub;
                    // todo: p_stage_default
                }
                int nx, nu, np_stage, ng_eq, ng_ineq;
                EvalBF L;
                EvalBF RSQrq;
                EvalBF rq;
                EvalBF g_equality;
                EvalBF g_inequality;
                EvalBF Ggt_equality;
                EvalBF Ggt_inequality;
                EvalBF BAbt;
                EvalBF b;
                vd lb;
                vd ub;
                vd p_stage_default;
            };
            std::vector<std::shared_ptr<evaluation_quantities>>  evaluation_quantities_ptr_;
            vd p_global_default_;
            int np_global_;
            int horizon_length_;
            std::vector<vd> x_initial;
            std::vector<vd> u_initial;
        };
    } // namespace casadi
} // namespace fatrop