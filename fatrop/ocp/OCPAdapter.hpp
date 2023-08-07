/*
 * Fatrop - A fast trajectory optimization solver
 * Copyright (C) 2022, 2023 Lander Vanroye, KU Leuven. All rights reserved.
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
#ifndef OCPEVALUATORINCLUDED
#define OCPEVALUATORINCLUDED
#include "OCPKKT.hpp"
#include "OCPAbstract.hpp"
#include "solver/FatropData.hpp"
#include "OCP.hpp"
#include <memory>
#include "auxiliary/Common.hpp"
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP->name)
#define AUXMACRO(type, name, suffix) type name##suffix = ((type)OCP->aux.name)
#define SOLVERMACRO(type, name, suffix) type name##suffix = ((type)name)
#ifdef ENABLE_MULTITHREADING
#include <omp.h>
#endif

// Example: you can make a for loop parallel with the following code. But to be further investigated to do this in an efficient way...
// #ifdef ENABLE_MULTITHREADING
// #pragma omp parallel for
// <your for loop>
// #endif

namespace fatrop
{
    class BFGSUpdater
    {
    public:
        BFGSUpdater(const int m) : m(m), Bk_prev(m, m), tmp1(m, 2), tmp2(m, 2), vk(m), yk_tilde(m) {}
        int update(MAT *Bip1, VEC *si, VEC *yi)
        {
            MAT *Bk_prev_p = Bk_prev;
            VEC *vk_p = vk;
            VEC *yk_tilde_p = yk_tilde;

            // powell update formula
            // compute update
            double sts = DOT(m, si, 0, si, 0);
            double sty = DOT(m, si, 0, yi, 0);
            double yty = DOT(m, yi, 0, yi, 0);
            if (sts == 0.0 || first_time)
            {
                skips = 0;
                reset();
                GECP(m, m, Bk_prev_p, 0, 0, Bip1, 0, 0);
                return 0;
            }
            GEMV_N(m, m, 1.0, Bk_prev_p, 0, 0, si, 0, 0.0, vk_p, 0, vk_p, 0);
            double stv = DOT(m, si, 0, vk_p, 0);
            double beta = -1.0 / stv;
#define SKIPPING
#define ALT_RANK1
#ifdef SKIPPING
            double alpha_tilde = 1.0 / sty;
            AXPBY(m, 1.0, yi, 0, 0.0, vk_p, 0, yk_tilde_p, 0);
            if (sty <= 1e-8 * std::sqrt(yty) * std::sqrt(sts))
            {
                int ret = 0;
                if (skips >= 1)
                {
                    reset(1.0);
                    ret = 1;
                    // std::cout << "resetting Bk" << std::endl;
                }
                GECP(m, m, Bk_prev_p, 0, 0, Bip1, 0, 0);
                skips++;
                return ret;
            }
#else
            double theta = sty > 0.2 * stv ? 1.0 : (0.8 * stv) / (stv - sty);
            AXPBY(m, theta, yi, 0, 1.0 - theta, vk_p, 0, yk_tilde_p, 0);
            double sty_tilde = DOT(m, si, 0, yk_tilde_p, 0);
            double alpha_tilde = 1.0 / sty_tilde;
#endif
#ifdef ALT_RANK1

            COLIN(m, yk_tilde_p, 0, tmp1, 0, 0);
            VECSC(m, alpha_tilde, yk_tilde_p, 0);
            COLIN(m, yk_tilde_p, 0, tmp2, 0, 0);

            COLIN(m, vk_p, 0, tmp1, 0, 1);
            VECSC(m, beta, vk_p, 0);
            COLIN(m, vk_p, 0, tmp2, 0, 1);

            SYRK_LN(m, 2, 1.0, tmp1, 0, 0, tmp2, 0, 0, 1.0, Bk_prev_p, 0, 0, Bip1, 0, 0);
            TRTR_L(m, Bip1, 0, 0, Bip1, 0, 0);
#else
            GER(m, m, alpha_tilde, yk_tilde_p, 0, yk_tilde_p, 0, Bk_prev_p, 0, 0, Bip1, 0, 0);
            GER(m, m, beta, vk_p, 0, vk_p, 0, Bip1, 0, 0, Bip1, 0, 0);
#endif
            // save the previous Bk
            GECP(m, m, Bip1, 0, 0, Bk_prev_p, 0, 0);
            skips = 0;
            return 0;
        }
        void reset(double alpha = 1.0, bool reset_skips = true)
        {
            // std:: cout << "resetting Bk" << std::endl;
            MAT *Bk_prev_p = Bk_prev;
            // identity matrix for B0
            GESE(m, m, 0.0, Bk_prev_p, 0, 0);
            DIARE(m, alpha, Bk_prev_p, 0, 0);
        }
        const int m;
        MATBF Bk_prev;
        MATBF tmp1;
        MATBF tmp2;
        VECBF vk;
        VECBF yk_tilde;
        int skips = 0;
        bool first_time = true;
    };
    class OCPBFGSUpdater : public BFGSUpdater
    {
    public:
        OCPBFGSUpdater(int nu, int nx, int nxp1, int ng, int ng_ineq, bool first, bool last) : BFGSUpdater(nu + nx), nu(nu), nx(nx), nxp1(nxp1), ng(ng), ng_ineq(ng_ineq), first(first), last(last), BAt_prev(nu + nx, nxp1), Gt_prev(nu + nx, ng), Gt_ineq_prev(nu + nx, ng_ineq), ux_prev(nu + nx), grad_obj_prev(nu + nx), s(nu + nx), y(nu + nx) {reset();}
        int update(MAT *Bkp1, VEC *ux, int a_ux, VEC *grad_obj, int a_grad_obj, MAT *BAbt, VEC *lam_dyn, int a_lam_dyn, MAT *Ggt, VEC *lam_eq, int a_lam_eq, MAT *Ggt_ineq, VEC *lam_ineq, int a_lam_ineq)
        {
            VEC *ux_prev_p = ux_prev;
            VEC *grad_obj_prev_p = grad_obj_prev;
            VEC *s_p = s;
            VEC *y_p = y;
            MAT *BAt_prev_p = BAt_prev;
            MAT *Gt_prev_p = Gt_prev;
            MAT *Gt_ineq_prev_p = Gt_ineq_prev;

            // compute s
            AXPBY(nu + nx, 1.0, ux, a_ux, -1.0, ux_prev_p, 0, s_p, 0);

            // compute y
            AXPBY(nu + nx, 1.0, grad_obj, a_grad_obj, 0.0, y_p, 0, y_p, 0);
            if (!last)
            {
                // contribution from dynamics
                GEMV_N(nu + nx, nxp1, 1.0, BAbt, 0, 0, lam_dyn, a_lam_dyn, 1.0, y_p, 0, y_p, 0);
            }
            // contribution from equality constraints
            GEMV_N(nu + nx, ng, 1.0, Ggt, 0, 0, lam_eq, a_lam_eq, 1.0, y_p, 0, y_p, 0);
            // contribution from inequality constraints
            GEMV_N(nu + nx, ng_ineq, 1.0, Ggt_ineq, 0, 0, lam_ineq, a_lam_ineq, 1.0, y_p, 0, y_p, 0);
            // save in last row
            ROWIN(nu + nx, 1.0, y_p, 0, Bkp1, nu + nx, 0);
            if (!first_time)
            {
                AXPBY(nu + nx, -1.0, grad_obj_prev_p, 0, 1.0, y_p, 0, y_p, 0);
                if (!last)
                {
                    // contribution from dynamics
                    GEMV_N(nu + nx, nxp1, -1.0, BAt_prev_p, 0, 0, lam_dyn, a_lam_dyn, 1.0, y_p, 0, y_p, 0);
                }
                // contribution from equality constraints
                GEMV_N(nu + nx, ng, -1.0, Gt_prev_p, 0, 0, lam_eq, a_lam_eq, 1.0, y_p, 0, y_p, 0);
                // contribution from inequality constraints
                GEMV_N(nu + nx, ng_ineq, -1.0, Gt_ineq_prev_p, 0, 0, lam_ineq, a_lam_ineq, 1.0, y_p, 0, y_p, 0);
            }
            // call BFGS update
            int ret = BFGSUpdater::update(Bkp1, s_p, y_p);
            // save ux and grad_obj
            VECCP(nu + nx, ux, a_ux, ux_prev_p, 0);
            VECCP(nu + nx, grad_obj, a_grad_obj, grad_obj_prev_p, 0);
            // save BAt, Gt, Gt_ineq
            GECP(nu + nx, nxp1, BAbt, 0, 0, BAt_prev_p, 0, 0);
            GECP(nu + nx, ng, Ggt, 0, 0, Gt_prev_p, 0, 0);
            GECP(nu + nx, ng_ineq, Ggt_ineq, 0, 0, Gt_ineq_prev_p, 0, 0);
            first_time = false;
            // blasfeo_print_dmat(nu+nx+1, nu+nx, Bkp1, 0, 0);
            return ret;
        }
        void reset()
        {
            first_time = true;
            BFGSUpdater::reset();
        }
        const int nu, nx, nxp1, ng, ng_ineq;
        const bool first, last;
        MATBF BAt_prev;
        MATBF Gt_prev;
        MATBF Gt_ineq_prev;
        VECBF ux_prev;
        VECBF grad_obj_prev;
        VECBF s;
        VECBF y;
    };
    class OCPAdapter : public OCP // public OCP -> also include KKTmemory, OCPDims, ...
    {
    public:
        OCPAdapter(const std::shared_ptr<OCPAbstract> &ocptempl_, const std::shared_ptr<FatropOptions> &options) : K(ocptempl_->get_horizon_length()),
                                                                                                                   nuexpr(TransformRange<fatrop_int>(0, K, [&ocptempl_](fatrop_int k)
                                                                                                                                                     { return ocptempl_->get_nuk(k); })),
                                                                                                                   nxexpr(TransformRange<fatrop_int>(0, K, [&ocptempl_](fatrop_int k)
                                                                                                                                                     { return ocptempl_->get_nxk(k); })),
                                                                                                                   ngexpr(TransformRange<fatrop_int>(0, K, [&ocptempl_](fatrop_int k)
                                                                                                                                                     { return ocptempl_->get_ngk(k); })),
                                                                                                                   ngineqexpr(TransformRange<fatrop_int>(0, K, [&ocptempl_](fatrop_int k)
                                                                                                                                                         { return ocptempl_->get_ng_ineq_k(k); })),
                                                                                                                   nstageparamsexpr(TransformRange<fatrop_int>(0, K, [&ocptempl_](fatrop_int k)
                                                                                                                                                               { return ocptempl_->get_n_stage_params_k(k); })),
                                                                                                                   offs_stageparams(offsets(nstageparamsexpr)), stageparams(sum(nstageparamsexpr), 0.0), globalparams(ocptempl_->get_n_global_params(), 0.0), options(options), ocptempl(ocptempl_)
        {
            #ifdef ENABLE_MULTITHREADING
            // check if environment variable OMP_NUM_THREADS is set
            if (getenv("OMP_NUM_THREADS") == NULL)
            {
                throw std::runtime_error("Environment variable OMP_NUM_THREADS is not set. Please set it to the number of threads you want to use or compile fatrop with multithreading disabled.");
            } 
            #endif

            // initialize the default parameters
            ocptempl_->get_default_global_params(globalparams.data());
            fatrop_int offs = 0;
            for (fatrop_int k = 0; k < K; k++)
            {
                ocptempl_->get_default_stage_paramsk(stageparams.data() + offs, k);
                offs += ocptempl_->get_n_stage_params_k(k);
            }
            for (fatrop_int k = 0; k < K; k++)
                OCPBFGS_updaters.emplace_back(ocptempl_->get_nuk(k), ocptempl_->get_nxk(k), K == K - 1 ? 0 : ocptempl_->get_nxk(k + 1), ocptempl_->get_ngk(k), ocptempl_->get_ng_ineq_k(k), k == 0, k == K - 1);
            // initialize gradbuf
            for (fatrop_int k = 0; k < K; k++)
                gradbuf.emplace_back(ocptempl_->get_nuk(k) + ocptempl_->get_nxk(k));
            x_dummy = std::vector<double>(maxel(nxexpr), 0.0);
            options->register_option(BooleanOption("bfgs", "bfgs Hessian approximation", &bfgs, false));
        }
        void reset()
        {
            if (bfgs)
            {
                for (auto &updater : OCPBFGS_updaters)
                    updater.reset();
            }
        }
        fatrop_int eval_lag_hess(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &lam) override;
        fatrop_int eval_constr_jac(
            OCPKKTMemory *OCP,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars) override;
        fatrop_int eval_contr_viol(
            OCPKKTMemory *OCP,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            FatropVecBF &constraint_violation) override;
        fatrop_int eval_obj_grad(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            FatropVecBF &gradient) override;
        fatrop_int eval_obj(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            double &res);
        fatrop_int integrate_dynamics(
            OCPKKTMemory *OCP,
            const fatrop_int k,
            const FatropVecBF &uk,
            const FatropVecBF &xk,
            FatropVecBF &xkp1);
        OCPDims get_ocp_dims() const override
        {
            return OCPDims(ocptempl->get_horizon_length(), nuexpr, nxexpr, ngexpr, ngineqexpr, nstageparamsexpr, ocptempl->get_n_global_params());
        }

    public:
        void set_parameters(const std::vector<double> &stage_params_in, const std::vector<double> &global_params_in) override;
        void set_initial_sol_guess(const std::shared_ptr<FatropData> &fatropdata, std::vector<double> &initial_u, std::vector<double> &initial_x) override;
        void get_solution(const std::shared_ptr<FatropData> &fatropdata, std::vector<double> &u, std::vector<double> &x) override;
        double *get_global_parameters()
        {
            if (globalparams.size() == 0)
                return nullptr;
            return globalparams.data();
        }
        double *get_stage_parameters()
        {
            if (stageparams.size() == 0)
                return nullptr;
            return stageparams.data();
        }
        std::vector<double> &get_global_parameters_vec()
        {
            return globalparams;
        }
        std::vector<double> &get_stage_parameters_vec()
        {
            return stageparams;
        }
        fatrop_int get_bounds(
            FatropVecBF &lower,
            FatropVecBF &upper) const override
        {
            fatrop_int offs = 0;
            double *lower_p = ((VEC *)lower)->pa;
            double *upper_p = ((VEC *)upper)->pa;
            for (fatrop_int k = 0; k < K; k++)
            {
                ocptempl->get_boundsk(lower_p + offs, upper_p + offs, k);
                offs += ocptempl->get_ng_ineq_k(k);
            }
            return 0;
        };
        fatrop_int get_initial_sol_guess(
            FatropVecBF &initial) const override
        {
            fatrop_int offs = 0;
            for (fatrop_int k = 0; k < K - 1; k++)
            {
                ocptempl->get_initial_uk(((VEC *)initial)->pa + offs, k);
                offs += ocptempl->get_nuk(k);
                ocptempl->get_initial_xk(((VEC *)initial)->pa + offs, k);
                offs += ocptempl->get_nxk(k);
            }
            ocptempl->get_initial_xk(((VEC *)initial)->pa + offs, K - 1);
            return 0;
        }
        // virtual fatrop_int GetDefaultParams(
        //     FatropOptions &params)
        // {
        // };

    public:
        fatrop_int K;
        FatropVector<fatrop_int> nuexpr;
        FatropVector<fatrop_int> nxexpr;
        FatropVector<fatrop_int> ngexpr;
        FatropVector<fatrop_int> ngineqexpr;
        FatropVector<fatrop_int> nstageparamsexpr;
        FatropVector<fatrop_int> offs_stageparams;
        std::vector<double> stageparams;
        std::vector<double> globalparams;
        std::vector<double> x_dummy;
        std::vector<OCPBFGSUpdater> OCPBFGS_updaters;
        std::vector<VECBF> gradbuf;
        std::shared_ptr<FatropOptions> options;
        bool bfgs;

    private:
        std::shared_ptr<OCPAbstract> ocptempl;
    };
} // namespace fatrop

#endif // OCPEVALUATORINCLUDED