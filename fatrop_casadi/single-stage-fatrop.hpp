#pragma once
#include <casadi/casadi.hpp>
#include <vector>
#include <memory>
#include "ocp/StageOCPApplication.hpp"
#include "function-evaluation.hpp"
#include "single-stage-problem.hpp"
/**
 *  Wraps a SingleStageProblem into fatrop::StageOCP
 */
namespace fatrop
{
    namespace fatrop_casadi
    {
        // class SingleStageProblemFatrop : public StageOCP
        // {
        // public:
        //     SingleStageProblemFatrop(const SingleStageProblem &problem) : StageOCP(problem.dims.nu, problem.dims.nx, problem.initial.ng_equality, problem.middle.ng_equality, problem.terminal.ng_equality, problem.initial.ng_inequality, problem.middle.ng_inequality, problem.terminal.ng_inequality, problem.initial.np_stage, problem.dims.np_global, problem.dims.K){

        //                                                                   };
        //     // implementation of the OCPAbstract methods
        //     int eval_BAbtk(const double *states_kp1,
        //                    const double *inputs_k,
        //                    const double *states_k,
        //                    const double *stage_params_k,
        //                    const double *global_params,
        //                    MAT *res,
        //                    const int k) override
        //     {
        //         const double *arg[5];
        //         arg[0] = states_k;
        //         arg[1] = states_kp1;
        //         arg[2] = inputs_k;
        //         arg[3] = stage_params_k;
        //         arg[4] = global_params;
        //         BAbt(arg, res);
        //         return 0;
        //     };
        //     int eval_RSQrqtk(const double *objective_scale,
        //                      const double *inputs_k,
        //                      const double *states_k,
        //                      const double *lam_dyn_k,
        //                      const double *lam_eq_k,
        //                      const double *lam_ineq_k,
        //                      const double *stage_params_k,
        //                      const double *global_params,
        //                      MAT *res,
        //                      const int k) override
        //     {
        //         const double *arg[7];
        //         arg[0] = states_k;
        //         arg[1] = inputs_k;
        //         arg[2] = stage_params_k;
        //         arg[3] = global_params;
        //         arg[4] = lam_dyn_k;
        //         arg[5] = lam_eq_k;
        //         arg[6] = lam_ineq_k;
        //         if (k == 0)
        //             initial.RSQrq(arg, res);
        //         else if (k == K_ - 1)
        //             terminal.RSQrq(arg, res);
        //         else
        //             middle.RSQrq(arg, res);
        //         return 0;
        //     }
        //     int eval_Ggtk(
        //         const double *inputs_k,
        //         const double *states_k,
        //         const double *stage_params_k,
        //         const double *global_params,
        //         MAT *res,
        //         const int k) override
        //     {
        //         const double *arg[4];
        //         arg[0] = states_k;
        //         arg[1] = inputs_k;
        //         arg[2] = stage_params_k;
        //         arg[3] = global_params;
        //         if (k == 0)
        //             initial.Ggt_equality(arg, res);
        //         else if (k == K_ - 1)
        //             terminal.Ggt_equality(arg, res);
        //         else
        //             middle.Ggt_equality(arg, res);
        //         return 0;
        //     }
        //     int eval_Ggt_ineqk(
        //         const double *inputs_k,
        //         const double *states_k,
        //         const double *stage_params_k,
        //         const double *global_params,
        //         MAT *res,
        //         const int k)
        //     {
        //         const double *arg[4];
        //         arg[0] = states_k;
        //         arg[1] = inputs_k;
        //         arg[2] = stage_params_k;
        //         arg[3] = global_params;
        //         if (k == 0)
        //             initial.Ggt_inequality(arg, res);
        //         else if (k == K_ - 1)
        //             terminal.Ggt_inequality(arg, res);
        //         else
        //             middle.Ggt_inequality(arg, res);
        //         return 0;
        //     }
        //     int eval_bk(
        //         const double *states_kp1,
        //         const double *inputs_k,
        //         const double *states_k,
        //         const double *stage_params_k,
        //         const double *global_params,
        //         double *res,
        //         const int k) override
        //     {
        //         const double *arg[5];
        //         arg[0] = states_k;
        //         arg[1] = states_kp1;
        //         arg[2] = inputs_k;
        //         arg[3] = stage_params_k;
        //         arg[4] = global_params;
        //         b(arg, res);
        //         return 0;
        //     }
        //     int eval_gk(
        //         const double *inputs_k,
        //         const double *states_k,
        //         const double *stage_params_k,
        //         const double *global_params,
        //         double *res,
        //         const int k) override
        //     {
        //         const double *arg[4];
        //         arg[0] = states_k;
        //         arg[1] = inputs_k;
        //         arg[2] = stage_params_k;
        //         arg[3] = global_params;
        //         if (k == 0)
        //             initial.g_equality(arg, res);
        //         else if (k == K_ - 1)
        //             terminal.g_equality(arg, res);
        //         else
        //             middle.g_equality(arg, res);
        //         return 0;
        //     }
        //     int eval_gineqk(
        //         const double *inputs_k,
        //         const double *states_k,
        //         const double *stage_params_k,
        //         const double *global_params,
        //         double *res,
        //         const int k) override
        //     {
        //         const double *arg[4];
        //         arg[0] = states_k;
        //         arg[1] = inputs_k;
        //         arg[2] = stage_params_k;
        //         arg[3] = global_params;
        //         if (k == 0)
        //         initial.g_inequality(arg, res);
        //         else if (k == K_ - 1)
        //         terminal.g_inequality(arg, res);
        //         else
        //         middle.g_inequality(arg, res);
        //         return 0;
        //     }
        //     int eval_rqk(
        //         const double *objective_scale,
        //         const double *inputs_k,
        //         const double *states_k,
        //         const double *stage_params_k,
        //         const double *global_params,
        //         double *res,
        //         const int k) override
        //     {
        //         const double *arg[4];
        //         arg[0] = states_k;
        //         arg[1] = inputs_k;
        //         arg[2] = stage_params_k;
        //         arg[3] = global_params;
        //         if (k == 0)
        //         initial.rq(arg, res);
        //         else if (k == K_ - 1)
        //         terminal.rq(arg, res);
        //         else
        //         middle.rq(arg, res);
        //         return 0;
        //     }

        //     int eval_Lk(
        //         const double *objective_scale,
        //         const double *inputs_k,
        //         const double *states_k,
        //         const double *stage_params_k,
        //         const double *global_params,
        //         double *res,
        //         const int k) override
        //     {
        //         const double *arg[4];
        //         arg[0] = states_k;
        //         arg[1] = inputs_k;
        //         arg[2] = stage_params_k;
        //         arg[3] = global_params;
        //         if (k == 0)
        //         initial.L(arg, res);
        //         else if (k == K_ - 1)
        //         terminal.L(arg, res);
        //         else
        //         middle.L(arg, res);
        //         return 0;
        //     }
        //     int get_initial_xk(double *xk, const int k) const override
        //     {
        //         int nx = this->get_nxk(k);
        //         for (int i = 0; i < nx; i++)
        //         {
        //             xk[i] = x_initial[k][i];
        //         }
        //         return 0;
        //     };
        //     int get_initial_uk(double *uk, const int k) const override
        //     {
        //         int nu = this->get_nuk(k);
        //         for (int i = 0; i < nu; i++)
        //         {
        //             uk[i] = u_initial[k][i];
        //         }
        //         return 0;
        //     };
        //     int set_initial_xk(double *xk, const int k)
        //     {
        //         return 0;
        //     };
        //     int set_initial_uk(double *uk, const int k)
        //     {
        //         return 0;
        //     };
        //     int get_boundsk(double *lower, double *upper, const int k) const override
        //     {
        //         int ngineq = this->get_ng_ineq_k(k);
        //         const double *lower_ineq;
        //         const double *upper_ineq;
        //         if (k == 0)
        //         {
        //             lower_ineq = initial.lb.data();
        //             upper_ineq = initial.lb.data();
        //         }
        //         else if (k == K_ - 1)
        //         {
        //             lower_ineq = terminal.lb.data();
        //             upper_ineq = terminal.ub.data();
        //         }
        //         else
        //         {
        //             lower_ineq = middle.lb.data();
        //             upper_ineq = middle.ub.data();
        //         }

        //         for (int j = 0; j < ngineq; j++)
        //         {
        //             lower[j] = lower_ineq[j];
        //             upper[j] = upper_ineq[j];
        //         }
        //         return 0;
        //     };
        //     int get_default_stage_paramsk(double *stage_params_res, const int k) const override
        //     {
        //         return 0;
        //     }
        //     int get_default_global_params(double *global_params_res) const override
        //     {
        //         return 0;
        //     }

        // private:
        //     typedef std::vector<double> vd;
        //     struct evaluation_quantities
        //     {
        //         eval_bf L;
        //         eval_bf RSQrq;
        //         eval_bf rq;
        //         eval_bf g_equality;
        //         eval_bf g_inequality;
        //         eval_bf Ggt_equality;
        //         eval_bf Ggt_inequality;
        //         vd lb;
        //         vd ub;
        //     };
        //     eval_bf BAbt;
        //     eval_bf b;
        //     evaluation_quantities initial;
        //     evaluation_quantities middle;
        //     evaluation_quantities terminal;
        //     std::vector<vd> x_initial;
        //     std::vector<vd> u_initial;
        // };
    }
};