#pragma once
extern "C"
{
#include <blasfeo.h>
}
#include "fatrop_ustage_eval_abstract.hpp"
#include "casadi_fe.hpp"
namespace fatrop
{
    namespace spectool
    {
        class FatropuStageEvalCasadi : public FatropuStageEvalAbstract
        {
            public:
            FatropuStageEvalCasadi(const uStageQuantities &sq, const cs::Dict &opts)
            {
                auto optss = opts;
                bool expand = cs::get_from_dict(optss, "post_expand", true);
                bool jit = cs::get_from_dict(optss, "jit", false);
                // we handle jit ourself
                optss["jit"] = false;
                cs::Dict jit_options_ = cs::get_from_dict(optss, "jit_options", casadi::Dict({{"flags", "-Ofast -march=native -ffast-math"}}));
                K_ = sq.K;
                nu_ = sq.nu();
                nx_ = sq.nx();
                np_stage_ = sq.np_stage();
                np_global_ = sq.np_global();
                ng_eq_ = sq.ng_eq();
                ng_ineq_ = sq.ng_ineq();
                nxp1_ = sq.nxp1();
                cs::MX lam_g_equality = cs::MX::sym("lam_g_equality", ng_eq_);
                cs::MX lam_g_inequality = cs::MX::sym("lam_g_inequality", ng_ineq_);
                cs::MX lam_dynamics = cs::MX::sym("lam_dynamics", nxp1_);
                std::vector<cs::MX> ustage_syms{sq.u, sq.x, sq.p_stage, sq.p_global};
                auto ux = cs::MX::veccat({sq.u, sq.x});
                L_ = CasadiFEWrap(cs::Function("L", ustage_syms, {sq.L}, optss), expand, jit, jit_options_);
                auto rq_sym = cs::MX();
                auto RSQ_sym = cs::MX::hessian(sq.L, ux, rq_sym);
                rq_ = CasadiFEWrap(cs::Function("rq", ustage_syms, {cs::MX::densify(rq_sym)}, optss), expand, jit, jit_options_);
                // if (sq.nxp1 > 0)
                {
                    auto xp1 = cs::MX::sym("xp1", nxp1_);
                    auto b = -xp1 + sq.x_next;
                    auto BAbt = cs::MX::horzcat({cs::MX::jacobian(sq.x_next, ux), b}).T();
                    b_ = CasadiFEWrap(cs::Function("b", {sq.x, xp1, sq.u, sq.p_stage, sq.p_global},
                                                   {densify(b)}, optss),
                                      expand, jit, jit_options_);
                    BAbt_ = CasadiFEWrap(cs::Function("BAbt", {sq.x, xp1, sq.u, sq.p_stage, sq.p_global},
                                                      {cs::MX::densify(BAbt)}, optss),
                                         expand, jit, jit_options_);
                    auto rq_dyn_sym = cs::MX();
                    auto RSQ_dyn_sym = cs::MX::hessian(cs::MX::dot(sq.x_next, lam_dynamics), ux, rq_dyn_sym);
                    rq_sym += rq_dyn_sym;
                    RSQ_sym += RSQ_dyn_sym;
                }
                // if(sq.ng_eq > 0)
                {
                    auto g_eq = sq.g;
                    auto Ggt_equality = cs::MX::horzcat({cs::MX::jacobian(g_eq, ux), g_eq}).T();
                    Ggt_equality_ = CasadiFEWrap(cs::Function("Ggt_equality", ustage_syms, {cs::MX::densify(Ggt_equality)}, optss), expand, jit, jit_options_);
                    g_equality_ = CasadiFEWrap(cs::Function("g_equality", ustage_syms, {cs::MX::densify(g_eq)}, optss), expand, jit, jit_options_);
                    auto rq_eq = cs::MX();
                    auto RSQ_eq = cs::MX::hessian(cs::MX::dot(g_eq, lam_g_equality), ux, rq_eq);
                    rq_sym += rq_eq;
                    RSQ_sym += RSQ_eq;
                }
                // if(sq.ng_ineq>0)
                {
                    auto g_ineq = sq.g_ineq;
                    auto Ggt_inequality = cs::MX::horzcat({cs::MX::jacobian(g_ineq, ux), g_ineq}).T();
                    Ggt_inequality_ = CasadiFEWrap(cs::Function("Ggt_inequality", ustage_syms, {cs::MX::densify(Ggt_inequality)}, optss), expand, jit, jit_options_);
                    g_inequality_ = CasadiFEWrap(cs::Function("g_inequality", ustage_syms, {cs::MX::densify(g_ineq)}, optss), expand, jit, jit_options_);
                    auto rq_ineq = cs::MX();
                    auto RSQ_ineq = cs::MX::hessian(cs::MX::dot(g_ineq, lam_g_inequality), ux, rq_ineq);
                    rq_sym += rq_ineq;
                    RSQ_sym += RSQ_ineq;
                }
                RSQrqt_ = CasadiFEWrap(cs::Function("RSQrqt", {sq.u, sq.x, lam_dynamics, lam_g_equality, lam_g_inequality, sq.p_stage, sq.p_global}, {cs::MX::densify(cs::MX::horzcat({RSQ_sym, rq_sym}).T())}, optss), expand, jit, jit_options_);
                Ub_ = sq.ub.nonzeros();
                Lb_ = sq.lb.nonzeros();
            }

        public:
            virtual int get_nx() const { return nx_; };
            virtual int get_nu() const { return nu_; };
            virtual int get_ng() const { return ng_eq_; };
            virtual int get_n_stage_params() const { return np_stage_; };
            virtual int get_ng_ineq() const { return ng_ineq_; };
            virtual int eval_BAbt(
                const double *states_kp1,
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res)
            {
                const double *arg[5];
                arg[0] = states_k;
                arg[1] = states_kp1;
                arg[2] = inputs_k;
                arg[3] = stage_params_k;
                arg[4] = global_params;
                BAbt_.eval(arg, res);
                return 0;
            };
            virtual int eval_RSQrqt(
                const double *objective_scale,
                const double *inputs_k,
                const double *states_k,
                const double *lam_dyn_k,
                const double *lam_eq_k,
                const double *lam_eq_ineq_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res)
            {
                const double *arg[7];
                arg[0] = inputs_k;
                arg[1] = states_k;
                arg[2] = lam_dyn_k;
                arg[3] = lam_eq_k;
                arg[4] = lam_eq_ineq_k;
                arg[5] = stage_params_k;
                arg[6] = global_params;
                RSQrqt_.eval(arg, res);
                return 0;
            };
            virtual int eval_Ggt(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res)
            {
                const double *arg[4];
                arg[0] = inputs_k;
                arg[1] = states_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                Ggt_equality_.eval(arg, res);
                return 0;
            };
            virtual int eval_Ggt_ineq(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                MAT *res)
            {
                const double *arg[4];
                arg[0] = inputs_k;
                arg[1] = states_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                Ggt_inequality_.eval(arg, res);
                return 0;
            };
            virtual int eval_b(
                const double *states_kp1,
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res)
            {
                const double *arg[5];
                arg[0] = states_k;
                arg[1] = states_kp1;
                arg[2] = inputs_k;
                arg[3] = stage_params_k;
                arg[4] = global_params;
                b_.eval(arg, res);
                return 0;
            }
            virtual int eval_g(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res) 
                {

                const double *arg[4];
                arg[0] = inputs_k;
                arg[1] = states_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                g_equality_.eval(arg, res);
                return 0;
                }
            virtual int eval_gineq(
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res) 
                {
                const double *arg[4];
                arg[0] = inputs_k;
                arg[1] = states_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                g_inequality_.eval(arg, res);
                return 0;

                }
            virtual int eval_rq(
                const double *objective_scale,
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res) 
                {
                const double *arg[4];
                arg[0] = inputs_k;
                arg[1] = states_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                rq_.eval(arg, res);
                return 0;

                }
            virtual int eval_L(
                const double *objective_scale,
                const double *inputs_k,
                const double *states_k,
                const double *stage_params_k,
                const double *global_params,
                double *res) 
                {
                const double *arg[4];
                arg[0] = inputs_k;
                arg[1] = states_k;
                arg[2] = stage_params_k;
                arg[3] = global_params;
                L_.eval(arg, res);
                return 0;
                }
            virtual int get_bounds(double *lower, double *upper) const 
            {
                for (int i = 0; i < ng_ineq_; i++)
                {
                    lower[i] = Lb_[i];
                    upper[i] = Ub_[i];
                }
                return 0;

            }
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