// Header-only runtime for the proxy -> Nlp<OcpType> codegen (see ocp_codegen.py).
//
// A generated <name>_ocp.cpp wires CasADi-generated C functions (one set per
// node "template") into the data structures here, and this header provides:
//   * CFun           : a CasADi C function + reusable arg/res buffers
//   * TemplateData   : per-template metadata + the CFuns for each eval kind
//   * GeneratedNlp   : a generic fatrop::Nlp<OcpType> that dispatches node k to
//                      node_to_template[k] and packs results into blasfeo
//   * Solver         : cold-start driver (fresh IpAlgBuilder per solve)
//
// Fatrop conventions (matching the CasADi fatrop interface byte-for-byte):
//   * per-node decision vector z = [u; x]  (controls first)
//   * eval_* matrices are (nu+nx+1) x ncols; only the top (nu+nx) rows are
//     filled by the user (Fatrop recomputes the bottom "internal" row from
//     the constraint value vectors, exactly as the reference hand-written
//     example leaves it zero).
//   * dynamics defect b = -x_{k+1} + f(x,u)
//
// CasADi function I/O conventions the codegen must honor (inputs in this order):
//   b,g,gineq             : (x, u, param0, ...) -> single dense vector
//   jac    (multi-output) : (x, u, param0, ...) -> (babt?, ggt?, ggt_ineq?)
//   obj    (multi-output) : (x, u, param0, ...) -> (L, rq)
//   rsqrqt                : (x, u, lam_dyn, lam_eq, lam_ineq, obj_scale, params...)
//   x0,u0                 : (param0, param1, ...)
// The constraint *value* functions (b, g, gineq) are emitted separately so
// CasADi's .expand() doesn't trip over multi-output shared CallSX nodes from
// the big-function kernel (this happens with nvblox + ca.external inside).
// The *Jacobian* function (jac) groups every stagewise Jacobian block in one
// multi-output call — that matches the granularity of fatrop::Nlp::eval_constr_jac
// and is the main shape change vs. the prior OcpAbstract-based codegen.
// "?" outputs in jac are omitted when their dimension would be zero
// (terminal nodes skip ``babt``; templates without (in)equalities skip
// ``ggt``/``ggt_ineq``).
// Outputs (each is a single dense column-major buffer):
//   b        -> f                            (nx_next)
//   g        -> g_eq                         (ng_eq)
//   gineq    -> g_ineq                       (ng_ineq)
//   babt     -> (d f / d z)^T                (nu+nx, nx_next)
//   ggt      -> (d g_eq / d z)^T             (nu+nx, ng_eq)
//   ggt_ineq -> (d g_ineq / d z)^T           (nu+nx, ng_ineq)
//   L        -> stage cost                   (1)            [obj_scale applied in C++]
//   rq       -> d L / d z                    (nu+nx)        [obj_scale applied in C++]
//   rsqrqt   -> Lagrangian Hessian wrt z     (nu+nx, nu+nx) [obj_scale baked in]
//   x0/u0    -> initial guess                (nx)/(nu)
#pragma once

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <fatrop/fatrop.hpp>
#include "casadi_c_wrap.hpp"

namespace ocp_direct
{
    using fatrop::Index;
    using fatrop::Scalar;
    // NOTE: MAT/VEC are global macros (blasfeo_dmat/dvec), not fatrop types,
    // so they are used unqualified (do not `using fatrop::MAT`).
    using casadi_c_wrap::CasadiCFunctions;
    using casadi_c_wrap::CasadiCFunctionWrap;

    // One CasADi C function with reusable I/O buffers.
    struct CFun
    {
        CasadiCFunctionWrap wrap;
        std::vector<std::vector<double>> arg;
        std::vector<std::vector<double>> res;
        explicit CFun(const CasadiCFunctions &f) : wrap(f) { wrap.allocate_arg_res(arg, res); }
        void call() { wrap(arg, res); }
    };

    // Kinds the generated <name>_ocp.cpp wires up.  B / G / GINEQ are single-
    // output value Functions.  JAC is a multi-output Function with up to three
    // outputs (babt, ggt, ggt_ineq) in a fixed order, skipping any whose
    // dimension would be zero.  OBJ always emits two outputs (L, rq).
    enum Kind { B, G, GINEQ, JAC, OBJ, RSQRQT, X0, U0, NKIND };

    // Per-template metadata + (optional) CFun for each eval kind.
    struct TemplateData
    {
        int nx = 0, nu = 0, nx_next = 0, ng_eq = 0, ng_ineq = 0;
        bool has_dynamics = false;
        std::vector<double> lb, ub;          // size ng_ineq
        std::unique_ptr<CFun> fn[NKIND];     // null if absent

        // Which output index of JAC each Jacobian block lives at.  -1 = absent.
        // Set in the generated cpp.
        int jac_idx_babt = -1, jac_idx_ggt = -1, jac_idx_ggt_ineq = -1;
    };

    // Build ProblemDims<OcpType> from per-stage nx/nu/ng/ng_ineq.
    inline fatrop::ProblemDims<fatrop::OcpType> make_ocp_dims(
        const std::vector<TemplateData> &templates,
        const std::vector<int> &node_to_template)
    {
        const Index K = static_cast<Index>(node_to_template.size());
        std::vector<Index> nu(K), nx(K), ng(K), ng_ineq(K);
        for (Index k = 0; k < K; ++k)
        {
            const TemplateData &t = templates[node_to_template[k]];
            nu[k] = t.nu;
            nx[k] = t.nx;
            ng[k] = t.ng_eq;
            ng_ineq[k] = t.ng_ineq;
        }
        return fatrop::ProblemDims<fatrop::OcpType>(K, std::move(nu), std::move(nx),
                                                    std::move(ng), std::move(ng_ineq));
    }

    // Build NlpDims (number_of_variables / eq / ineq) from ProblemDims.
    inline fatrop::NlpDims make_nlp_dims(const fatrop::ProblemDims<fatrop::OcpType> &d)
    {
        Index nvar = 0, neq = 0, nineq = 0;
        for (Index k = 0; k < d.K; ++k)
        {
            nvar += d.number_of_controls[k] + d.number_of_states[k];
            neq += d.number_of_eq_constraints[k] + d.number_of_ineq_constraints[k];
            if (k != d.K - 1) neq += d.number_of_states[k + 1];
            nineq += d.number_of_ineq_constraints[k];
        }
        return fatrop::NlpDims(nvar, neq, nineq);
    }

    class GeneratedNlp : public fatrop::Nlp<fatrop::OcpType>
    {
    public:
        // Built by the generated .cpp.
        std::vector<TemplateData> templates;
        std::vector<int> node_to_template;
        std::vector<std::pair<int, int>> param_dims;   // (rows, cols) per param

        // Runtime state.
        std::vector<std::vector<double>> params;       // flat (rows*cols) per param
        std::vector<std::vector<double>> x0_cache;     // per node (nx)
        std::vector<std::vector<double>> u0_cache;     // per node (nu)
        bool init_overridden = false;

        // Dim caches (filled in `finalize()` once `templates`/`node_to_template`
        // are populated by the generated cpp).
        std::unique_ptr<fatrop::ProblemDims<fatrop::OcpType>> ocp_dims_;
        std::unique_ptr<fatrop::NlpDims> nlp_dims_;

        TemplateData &tmpl(Index k) { return templates[node_to_template[k]]; }
        const TemplateData &tmpl(Index k) const { return templates[node_to_template[k]]; }
        Index horizon() const { return static_cast<Index>(node_to_template.size()); }

        const fatrop::NlpDims &nlp_dims() const override { return *nlp_dims_; }
        const fatrop::ProblemDims<fatrop::OcpType> &problem_dims() const override
        {
            return *ocp_dims_;
        }

        // --- helpers ---
        static void copy_in(std::vector<double> &dst, const Scalar *src)
        {
            if (src) std::copy(src, src + dst.size(), dst.begin());
            else std::fill(dst.begin(), dst.end(), 0.0);
        }
        void fill_xu_p(CFun &cf, const Scalar *x, const Scalar *u)
        {
            copy_in(cf.arg[0], x);
            copy_in(cf.arg[1], u);
            for (size_t i = 0; i < params.size(); ++i) cf.arg[2 + i] = params[i];
        }
        void fill_p_only(CFun &cf)
        {
            for (size_t i = 0; i < params.size(); ++i) cf.arg[i] = params[i];
        }
        static void pack(MAT *res, int rows, int cols, const std::vector<double> &buf)
        {
            fatrop::blasfeo_gese_wrap(res->m, res->n, 0.0, res, 0, 0);
            if (rows > 0 && cols > 0)
                fatrop::blasfeo_pack_mat_wrap(rows, cols, const_cast<double *>(buf.data()),
                                              rows, res, 0, 0);
        }

        // Build the dim caches.  Call after `templates` / `node_to_template`
        // are populated.
        void finalize_dims()
        {
            ocp_dims_ = std::make_unique<fatrop::ProblemDims<fatrop::OcpType>>(
                make_ocp_dims(templates, node_to_template));
            nlp_dims_ = std::make_unique<fatrop::NlpDims>(make_nlp_dims(*ocp_dims_));
        }

        // ------------------------------------------------------------ Nlp API

        Index eval_lag_hess(const fatrop::ProblemInfo<fatrop::OcpType> &info,
                            const Scalar objective_scale, const fatrop::VecRealView &primal_x,
                            const fatrop::VecRealView &primal_s,
                            const fatrop::VecRealView &lam,
                            fatrop::Hessian<fatrop::OcpType> &hess) override
        {
            const Scalar *primal_x_p = primal_x.data();
            const Scalar *lam_p = lam.data();
            for (Index k = 0; k < info.dims.K; ++k)
            {
                TemplateData &t = tmpl(k);
                CFun &cf = *t.fn[RSQRQT];
                const Scalar *u_k = primal_x_p + info.offsets_primal_u[k];
                const Scalar *x_k = primal_x_p + info.offsets_primal_x[k];
                const Scalar *lam_dyn_k =
                    (k != info.dims.K - 1) ? lam_p + info.offsets_g_eq_dyn[k] : nullptr;
                const Scalar *lam_eq_k = lam_p + info.offsets_g_eq_path[k];
                const Scalar *lam_ineq_k = lam_p + info.offsets_g_eq_slack[k];

                // CasADi sig: (x, u, lam_dyn, lam_eq, lam_ineq, obj_scale, params...)
                copy_in(cf.arg[0], x_k);
                copy_in(cf.arg[1], u_k);
                copy_in(cf.arg[2], lam_dyn_k);
                copy_in(cf.arg[3], lam_eq_k);
                copy_in(cf.arg[4], lam_ineq_k);
                cf.arg[5][0] = objective_scale;
                for (size_t i = 0; i < params.size(); ++i) cf.arg[6 + i] = params[i];
                cf.call();
                pack(&hess.RSQrqt[k].mat(), t.nu + t.nx, t.nu + t.nx, cf.res[0]);
            }
            return 0;
        }

        Index eval_constr_jac(const fatrop::ProblemInfo<fatrop::OcpType> &info,
                              const fatrop::VecRealView &primal_x,
                              const fatrop::VecRealView &primal_s,
                              fatrop::Jacobian<fatrop::OcpType> &jac) override
        {
            const Scalar *primal_x_p = primal_x.data();
            for (Index k = 0; k < info.dims.K; ++k)
            {
                TemplateData &t = tmpl(k);
                if (!t.fn[JAC]) continue;       // template has no Jacobian blocks
                CFun &cf = *t.fn[JAC];
                fill_xu_p(cf, primal_x_p + info.offsets_primal_x[k],
                          primal_x_p + info.offsets_primal_u[k]);
                cf.call();
                if (t.jac_idx_babt >= 0)
                    pack(&jac.BAbt[k].mat(), t.nu + t.nx, t.nx_next,
                         cf.res[t.jac_idx_babt]);
                if (t.jac_idx_ggt >= 0)
                    pack(&jac.Gg_eqt[k].mat(), t.nu + t.nx, t.ng_eq,
                         cf.res[t.jac_idx_ggt]);
                if (t.jac_idx_ggt_ineq >= 0)
                    pack(&jac.Gg_ineqt[k].mat(), t.nu + t.nx, t.ng_ineq,
                         cf.res[t.jac_idx_ggt_ineq]);
            }
            return 0;
        }

        Index eval_constraint_violation(const fatrop::ProblemInfo<fatrop::OcpType> &info,
                                        const fatrop::VecRealView &primal_x,
                                        const fatrop::VecRealView &primal_s,
                                        fatrop::VecRealView &res) override
        {
            const Scalar *primal_x_p = primal_x.data();
            Scalar *res_p = res.data();
            for (Index k = 0; k < info.dims.K; ++k)
            {
                TemplateData &t = tmpl(k);
                const Scalar *x_k = primal_x_p + info.offsets_primal_x[k];
                const Scalar *u_k = primal_x_p + info.offsets_primal_u[k];
                if (t.fn[B])
                {
                    CFun &cf = *t.fn[B];
                    fill_xu_p(cf, x_k, u_k);
                    cf.call();
                    Scalar *b_out = res_p + info.offsets_g_eq_dyn[k];
                    const Scalar *x_kp1 = primal_x_p + info.offsets_primal_x[k + 1];
                    for (int i = 0; i < t.nx_next; ++i) b_out[i] = -x_kp1[i] + cf.res[0][i];
                }
                if (t.fn[G])
                {
                    CFun &cf = *t.fn[G];
                    fill_xu_p(cf, x_k, u_k);
                    cf.call();
                    std::copy(cf.res[0].begin(), cf.res[0].end(),
                              res_p + info.offsets_g_eq_path[k]);
                }
                if (t.fn[GINEQ])
                {
                    CFun &cf = *t.fn[GINEQ];
                    fill_xu_p(cf, x_k, u_k);
                    cf.call();
                    std::copy(cf.res[0].begin(), cf.res[0].end(),
                              res_p + info.offsets_g_eq_slack[k]);
                }
            }
            // add -s to the slack constraints (g_ineq - s = 0).
            res.block(info.number_of_g_eq_slack, info.offset_g_eq_slack) =
                res.block(info.number_of_g_eq_slack, info.offset_g_eq_slack) -
                primal_s.block(info.number_of_g_eq_slack, 0);
            return 0;
        }

        Index eval_objective_gradient(const fatrop::ProblemInfo<fatrop::OcpType> &info,
                                      const Scalar objective_scale,
                                      const fatrop::VecRealView &primal_x,
                                      const fatrop::VecRealView &primal_s,
                                      fatrop::VecRealView &grad_x,
                                      fatrop::VecRealView &grad_s) override
        {
            grad_s.block(info.number_of_g_eq_slack, 0) = 0;
            const Scalar *primal_x_p = primal_x.data();
            Scalar *grad_x_p = grad_x.data();
            for (Index k = 0; k < info.dims.K; ++k)
            {
                TemplateData &t = tmpl(k);
                CFun &cf = *t.fn[OBJ];
                fill_xu_p(cf, primal_x_p + info.offsets_primal_x[k],
                          primal_x_p + info.offsets_primal_u[k]);
                cf.call();
                // out[1] = rq, length nu+nx, ordered [u; x] (matches z = [u; x]).
                Scalar *out = grad_x_p + info.offsets_primal_u[k];
                const auto &rq = cf.res[1];
                for (int i = 0; i < t.nu + t.nx; ++i) out[i] = objective_scale * rq[i];
            }
            return 0;
        }

        Index eval_objective(const fatrop::ProblemInfo<fatrop::OcpType> &info,
                             const Scalar objective_scale,
                             const fatrop::VecRealView &primal_x,
                             const fatrop::VecRealView &primal_s, Scalar &res) override
        {
            res = 0;
            const Scalar *primal_x_p = primal_x.data();
            for (Index k = 0; k < info.dims.K; ++k)
            {
                CFun &cf = *tmpl(k).fn[OBJ];
                fill_xu_p(cf, primal_x_p + info.offsets_primal_x[k],
                          primal_x_p + info.offsets_primal_u[k]);
                cf.call();
                res += objective_scale * cf.res[0][0];
            }
            return 0;
        }

        Index get_bounds(const fatrop::ProblemInfo<fatrop::OcpType> &info,
                         fatrop::VecRealView &lower_bounds,
                         fatrop::VecRealView &upper_bounds) override
        {
            if (info.number_of_slack_variables == 0) return 0;
            Scalar *lb_p = lower_bounds.data();
            Scalar *ub_p = upper_bounds.data();
            for (Index k = 0; k < info.dims.K; ++k)
            {
                const TemplateData &t = tmpl(k);
                Scalar *lb_k = lb_p + info.offsets_slack[k];
                Scalar *ub_k = ub_p + info.offsets_slack[k];
                for (int i = 0; i < t.ng_ineq; ++i)
                {
                    lb_k[i] = t.lb[i];
                    ub_k[i] = t.ub[i];
                }
            }
            return 0;
        }

        Index get_initial_primal(const fatrop::ProblemInfo<fatrop::OcpType> &info,
                                 fatrop::VecRealView &primal_x) override
        {
            Scalar *primal_x_p = primal_x.data();
            for (Index k = 0; k < info.dims.K; ++k)
            {
                const TemplateData &t = tmpl(k);
                const auto &xv = x0_cache[k];
                const auto &uv = u0_cache[k];
                std::copy(uv.begin(), uv.end(), primal_x_p + info.offsets_primal_u[k]);
                std::copy(xv.begin(), xv.end(), primal_x_p + info.offsets_primal_x[k]);
                (void)t;
            }
            return 0;
        }

        void get_primal_damping(const fatrop::ProblemInfo<fatrop::OcpType> &info,
                                fatrop::VecRealView &damping) override
        {
            damping = 0;
        }

        void apply_jacobian_s_transpose(const fatrop::ProblemInfo<fatrop::OcpType> &info,
                                        const fatrop::VecRealView &multipliers,
                                        const Scalar alpha,
                                        const fatrop::VecRealView &y,
                                        fatrop::VecRealView &out) override
        {
            out = alpha * y;
            out.block(info.number_of_slack_variables, 0) =
                out.block(info.number_of_slack_variables, 0) -
                multipliers.block(info.number_of_slack_variables, info.offset_g_eq_slack);
        }

        // Size the x0/u0 caches once (call before the first solve). Idempotent.
        void allocate_caches()
        {
            Index K = horizon();
            x0_cache.resize(K);
            u0_cache.resize(K);
            for (Index k = 0; k < K; ++k)
            {
                TemplateData &t = tmpl(k);
                if ((int)x0_cache[k].size() != t.nx) x0_cache[k].assign(t.nx, 0.0);
                if ((int)u0_cache[k].size() != t.nu) u0_cache[k].assign(t.nu, 0.0);
            }
        }

        // Refill x0/u0 from the x0/u0 CFuns and current params, unless an explicit
        // initial guess was set via ocp_set_initial. No (re)allocation: the caches
        // are pre-sized by allocate_caches() and only overwritten here.
        void fill_initial()
        {
            if (init_overridden) return;
            Index K = horizon();
            for (Index k = 0; k < K; ++k)
            {
                TemplateData &t = tmpl(k);
                std::fill(x0_cache[k].begin(), x0_cache[k].end(), 0.0);
                std::fill(u0_cache[k].begin(), u0_cache[k].end(), 0.0);
                if (t.fn[X0])
                {
                    fill_p_only(*t.fn[X0]);
                    t.fn[X0]->call();
                    std::copy(t.fn[X0]->res[0].begin(), t.fn[X0]->res[0].end(),
                              x0_cache[k].begin());
                }
                if (t.fn[U0])
                {
                    fill_p_only(*t.fn[U0]);
                    t.fn[U0]->call();
                    std::copy(t.fn[U0]->res[0].begin(), t.fn[U0]->res[0].end(),
                              u0_cache[k].begin());
                }
            }
        }
    };

    // Driver. The interior-point algorithm is built ONCE (build() is where all
    // allocation happens) and reused: Fatrop's optimize() internally calls
    // reset() + initializer_->initialize(), so each call is a fresh cold solve
    // that re-reads get_initial_primal and the (mutated) parameters. solve() does
    // no dynamic allocation of its own — params, initial-guess caches and the
    // x_sol/u_sol read-back buffers are all pre-sized during setup().
    struct Solver
    {
        std::shared_ptr<GeneratedNlp> nlp;
        std::map<std::string, double> opt_d;
        std::map<std::string, long long> opt_i;
        std::map<std::string, bool> opt_b;

        // Built once in setup().
        fatrop::OptionRegistry options;
        std::unique_ptr<fatrop::IpAlgBuilder<fatrop::OcpType>> builder;
        std::shared_ptr<fatrop::IpAlgorithm<fatrop::OcpType>> ipalg;
        std::shared_ptr<fatrop::IpData<fatrop::OcpType>> ipdata;
        bool is_setup = false;

        // last-solve results (pre-sized read-back buffers)
        int iter_count = -1;
        double solve_time = 0.0;
        int ret_flag = -1;
        std::vector<std::vector<double>> x_sol;   // per node (nx)
        std::vector<std::vector<double>> u_sol;   // per node (nu)

        void setup()
        {
            nlp->finalize_dims();
            nlp->allocate_caches();
            builder = std::make_unique<fatrop::IpAlgBuilder<fatrop::OcpType>>(nlp);
            ipalg = builder->with_options_registry(&options).build();
            ipdata = builder->get_ipdata();
            // Options must be applied AFTER build() (the registry is populated then).
            for (auto &kv : opt_d) options.set_option<double>(kv.first, kv.second);
            for (auto &kv : opt_i)
                options.set_option<fatrop::Index>(kv.first, (fatrop::Index)kv.second);
            for (auto &kv : opt_b) options.set_option<bool>(kv.first, kv.second);
            Index K = nlp->horizon();
            x_sol.resize(K);
            u_sol.resize(K);
            for (Index k = 0; k < K; ++k)
            {
                x_sol[k].assign(nlp->tmpl(k).nx, 0.0);
                u_sol[k].assign(nlp->tmpl(k).nu, 0.0);
            }
            is_setup = true;
        }

        int solve()
        {
            if (!is_setup) setup();          // one-time
            nlp->fill_initial();             // no allocation (caches pre-sized)
            ret_flag = (int)ipalg->optimize();
            iter_count = (int)ipdata->iteration_number();
            solve_time = ipdata->timing_statistics().full_algorithm.elapsed();

            const auto &x = ipalg->solution_primal();
            const auto &info = ipalg->info();
            Index K = nlp->horizon();
            for (Index k = 0; k < K; ++k)
            {
                int nx = nlp->tmpl(k).nx;
                int nu = nlp->tmpl(k).nu;
                int ox = info.offsets_primal_x[k], ou = info.offsets_primal_u[k];
                for (int i = 0; i < nx; ++i) x_sol[k][i] = x(ox + i);
                for (int i = 0; i < nu; ++i) u_sol[k][i] = x(ou + i);
            }
            return ret_flag;
        }
    };
} // namespace ocp_direct
