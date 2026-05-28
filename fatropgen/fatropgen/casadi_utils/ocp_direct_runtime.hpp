// Header-only runtime for the proxy -> OcpAbstract codegen (see core/ocp_codegen.py).
//
// A generated <name>_ocp.cpp wires CasADi-generated C functions (one set per
// node "template") into the data structures here, and this header provides:
//   * CFun           : a CasADi C function + reusable arg/res buffers
//   * TemplateData   : per-template metadata + the CFuns for each eval kind
//   * GeneratedOcp   : a generic fatrop::OcpAbstract that dispatches node k to
//                      node_to_template[k] and packs results into blasfeo
//   * Solver         : cold-start driver (fresh IpAlgBuilder per solve)
//
// Fatrop conventions (matching the CasADi fatrop interface byte-for-byte):
//   * per-node decision vector z = [u; x]  (controls first)
//   * eval_* matrices are (nu+nx+1) x ncols; only the top (nu+nx) rows are
//     filled by the user (Fatrop recomputes the bottom "internal" row from
//     eval_b / eval_g / eval_gineq / eval_rq, exactly as the reference
//     hand-written example leaves it zero).
//   * dynamics defect b = -x_{k+1} + f(x,u)
//
// CasADi function I/O conventions the codegen must honor (inputs in this order):
//   b,babt,g,ggt,gineq,ggt_ineq,L,rq : (x, u, param0, param1, ...)
//   rsqrqt                            : (x, u, lam_dyn, lam_eq, lam_ineq, obj_scale, params...)
//   x0,u0                             : (param0, param1, ...)
// Outputs (single, dense, column-major):
//   b        -> f                        (nx_next)
//   babt     -> (d f / d z)^T            (nu+nx, nx_next)
//   g        -> g_eq                     (ng_eq)
//   ggt      -> (d g_eq / d z)^T         (nu+nx, ng_eq)
//   gineq    -> g_ineq                   (ng_ineq)
//   ggt_ineq -> (d g_ineq / d z)^T       (nu+nx, ng_ineq)
//   L        -> stage cost               (1)            [scaled by obj_scale in C]
//   rq       -> d L / d z                (nu+nx)        [scaled by obj_scale in C]
//   rsqrqt   -> Lagrangian Hessian wrt z (nu+nx, nu+nx) [obj_scale baked in]
//   x0/u0    -> initial guess            (nx)/(nu)
#pragma once

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <cassert>
#include <fatrop/fatrop.hpp>

// Construct a CFun from a CasADi-generated C function `name`. Expects the
// standard symbols emitted by `casadi::CodeGenerator` (name, name_n_in, ...).
#define MK_CFUN(name) std::make_unique<ocp_direct::CFun>( \
    &name, &name##_work, &name##_sparsity_in, &name##_sparsity_out, \
    &name##_n_in, &name##_n_out, &name##_incref, &name##_decref, \
    &name##_checkout, &name##_release)

namespace ocp_direct
{
    using fatrop::Index;
    using fatrop::Scalar;
    // NOTE: MAT/VEC are global macros (blasfeo_dmat/dvec), not fatrop types,
    // so they are used unqualified (do not `using fatrop::MAT`).

    // One CasADi C function with reusable I/O buffers. Owns the casadi memory
    // checkout for its lifetime (incref/checkout in ctor, release/decref in dtor).
    struct CFun
    {
        using casadi_int_t = long long int;
        using eval_t = int (*)(const double **, double **, casadi_int_t *, double *, int);
        using work_t = int (*)(casadi_int_t *, casadi_int_t *, casadi_int_t *, casadi_int_t *);
        using sparsity_t = const casadi_int_t *(*)(casadi_int_t);
        using count_t = casadi_int_t (*)();
        using void_t = void (*)();
        using checkout_t = int (*)();
        using release_t = void (*)(int);

        eval_t eval_;
        void_t decref_;
        release_t release_;
        int mem_id_;

        // User-facing buffers (one vector per input/output, dense storage).
        std::vector<std::vector<double>> arg;
        std::vector<std::vector<double>> res;
        // Scratch used by the casadi eval call.
        std::vector<const double *> arg_ptrs_;
        std::vector<double *> res_ptrs_;
        std::vector<casadi_int_t> iw_;
        std::vector<double> w_;

        CFun(eval_t eval, work_t work, sparsity_t sp_in, sparsity_t sp_out,
             count_t n_in, count_t n_out, void_t incref, void_t decref,
             checkout_t checkout, release_t release)
            : eval_(eval), decref_(decref), release_(release)
        {
            incref();
            mem_id_ = checkout();

            casadi_int_t ni = n_in(), no = n_out();
            arg.resize(ni);
            res.resize(no);
            for (casadi_int_t i = 0; i < ni; ++i)
            {
                const casadi_int_t *sp = sp_in(i);
                arg[i].assign(sp ? sp[0] * sp[1] : 0, 0.0);
            }
            for (casadi_int_t i = 0; i < no; ++i)
            {
                const casadi_int_t *sp = sp_out(i);
                res[i].assign(sp ? sp[0] * sp[1] : 0, 0.0);
            }

            casadi_int_t sz_arg = ni, sz_res = no, sz_iw = 0, sz_w = 0;
            work(&sz_arg, &sz_res, &sz_iw, &sz_w);
            arg_ptrs_.resize(sz_arg);
            res_ptrs_.resize(sz_res);
            iw_.resize(sz_iw);
            w_.resize(sz_w);
        }
        ~CFun()
        {
            release_(mem_id_);
            decref_();
        }
        CFun(const CFun &) = delete;
        CFun &operator=(const CFun &) = delete;

        void call()
        {
            for (size_t i = 0; i < arg.size(); ++i) arg_ptrs_[i] = arg[i].data();
            for (size_t i = 0; i < res.size(); ++i) res_ptrs_[i] = res[i].data();
            eval_(arg_ptrs_.data(), res_ptrs_.data(), iw_.data(), w_.data(), mem_id_);
        }
    };

    enum Kind { B, BABT, G, GGT, GINEQ, GGT_INEQ, L, RQ, RSQRQT, X0, U0, NKIND };

    // Per-template metadata + (optional) CFun for each eval kind.
    struct TemplateData
    {
        int nx = 0, nu = 0, nx_next = 0, ng_eq = 0, ng_ineq = 0;
        bool has_dynamics = false;
        std::vector<double> lb, ub;          // size ng_ineq
        std::unique_ptr<CFun> fn[NKIND];     // null if absent
    };

    class GeneratedOcp : public fatrop::OcpAbstract
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

        TemplateData &tmpl(Index k) { return templates[node_to_template[k]]; }
        const TemplateData &tmpl(Index k) const { return templates[node_to_template[k]]; }

        Index get_horizon_length() const override { return (Index)node_to_template.size(); }
        Index get_nx(const Index k) const override { return tmpl(k).nx; }
        Index get_nu(const Index k) const override { return tmpl(k).nu; }
        Index get_ng(const Index k) const override { return tmpl(k).ng_eq; }
        Index get_ng_ineq(const Index k) const override { return tmpl(k).ng_ineq; }

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

        // --- dynamics ---
        Index eval_b(const Scalar *x_kp1, const Scalar *u, const Scalar *x,
                     Scalar *res, const Index k) override
        {
            TemplateData &t = tmpl(k);
            CFun &cf = *t.fn[B];
            fill_xu_p(cf, x, u);
            cf.call();
            for (int i = 0; i < t.nx_next; ++i) res[i] = -x_kp1[i] + cf.res[0][i];
            return 0;
        }
        Index eval_BAbt(const Scalar *, const Scalar *u, const Scalar *x,
                        MAT *res, const Index k) override
        {
            TemplateData &t = tmpl(k);
            CFun &cf = *t.fn[BABT];
            fill_xu_p(cf, x, u);
            cf.call();
            pack(res, t.nu + t.nx, t.nx_next, cf.res[0]);
            return 0;
        }

        // --- equality constraints ---
        Index eval_g(const Scalar *u, const Scalar *x, Scalar *res, const Index k) override
        {
            TemplateData &t = tmpl(k);
            if (t.ng_eq == 0) return 0;
            CFun &cf = *t.fn[G];
            fill_xu_p(cf, x, u);
            cf.call();
            std::copy(cf.res[0].begin(), cf.res[0].end(), res);
            return 0;
        }
        Index eval_Ggt(const Scalar *u, const Scalar *x, MAT *res, const Index k) override
        {
            TemplateData &t = tmpl(k);
            if (t.ng_eq == 0) return 0;
            CFun &cf = *t.fn[GGT];
            fill_xu_p(cf, x, u);
            cf.call();
            pack(res, t.nu + t.nx, t.ng_eq, cf.res[0]);
            return 0;
        }

        // --- inequality constraints ---
        Index eval_gineq(const Scalar *u, const Scalar *x, Scalar *res, const Index k) override
        {
            TemplateData &t = tmpl(k);
            if (t.ng_ineq == 0) return 0;
            CFun &cf = *t.fn[GINEQ];
            fill_xu_p(cf, x, u);
            cf.call();
            std::copy(cf.res[0].begin(), cf.res[0].end(), res);
            return 0;
        }
        Index eval_Ggt_ineq(const Scalar *u, const Scalar *x, MAT *res, const Index k) override
        {
            TemplateData &t = tmpl(k);
            if (t.ng_ineq == 0) return 0;
            CFun &cf = *t.fn[GGT_INEQ];
            fill_xu_p(cf, x, u);
            cf.call();
            pack(res, t.nu + t.nx, t.ng_ineq, cf.res[0]);
            return 0;
        }
        Index get_bounds(Scalar *lower, Scalar *upper, const Index k) const override
        {
            const TemplateData &t = tmpl(k);
            for (int i = 0; i < t.ng_ineq; ++i) { lower[i] = t.lb[i]; upper[i] = t.ub[i]; }
            return 0;
        }

        // --- cost ---
        Index eval_L(const Scalar *obj_scale, const Scalar *u, const Scalar *x,
                     Scalar *res, const Index k) override
        {
            CFun &cf = *tmpl(k).fn[L];
            fill_xu_p(cf, x, u);
            cf.call();
            *res = obj_scale[0] * cf.res[0][0];
            return 0;
        }
        Index eval_rq(const Scalar *obj_scale, const Scalar *u, const Scalar *x,
                      Scalar *res, const Index k) override
        {
            TemplateData &t = tmpl(k);
            CFun &cf = *t.fn[RQ];
            fill_xu_p(cf, x, u);
            cf.call();
            for (int i = 0; i < t.nu + t.nx; ++i) res[i] = obj_scale[0] * cf.res[0][i];
            return 0;
        }

        // --- Lagrangian Hessian ---
        Index eval_RSQrqt(const Scalar *obj_scale, const Scalar *u, const Scalar *x,
                          const Scalar *lam_dyn, const Scalar *lam_eq,
                          const Scalar *lam_eq_ineq, MAT *res, const Index k) override
        {
            TemplateData &t = tmpl(k);
            CFun &cf = *t.fn[RSQRQT];
            // inputs: x, u, lam_dyn, lam_eq, lam_ineq, obj_scale, params...
            copy_in(cf.arg[0], x);
            copy_in(cf.arg[1], u);
            copy_in(cf.arg[2], lam_dyn);       // size nx_next (0 / null at terminal)
            copy_in(cf.arg[3], lam_eq);        // size ng_eq
            copy_in(cf.arg[4], lam_eq_ineq);   // size ng_ineq
            cf.arg[5][0] = obj_scale[0];
            for (size_t i = 0; i < params.size(); ++i) cf.arg[6 + i] = params[i];
            cf.call();
            pack(res, t.nu + t.nx, t.nu + t.nx, cf.res[0]);
            return 0;
        }

        // --- initial guess (const; reads cache prepared in prepare()) ---
        Index get_initial_xk(Scalar *xk, const Index k) const override
        {
            const auto &v = x0_cache[k];
            std::copy(v.begin(), v.end(), xk);
            return 0;
        }
        Index get_initial_uk(Scalar *uk, const Index k) const override
        {
            const auto &v = u0_cache[k];
            std::copy(v.begin(), v.end(), uk);
            return 0;
        }

        // Size the x0/u0 caches once (call before the first solve). Idempotent.
        void allocate_caches()
        {
            int K = get_horizon_length();
            x0_cache.resize(K);
            u0_cache.resize(K);
            for (int k = 0; k < K; ++k)
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
            int K = get_horizon_length();
            for (int k = 0; k < K; ++k)
            {
                TemplateData &t = tmpl(k);
                std::fill(x0_cache[k].begin(), x0_cache[k].end(), 0.0);
                std::fill(u0_cache[k].begin(), u0_cache[k].end(), 0.0);
                if (t.fn[X0]) { fill_p_only(*t.fn[X0]); t.fn[X0]->call();
                                std::copy(t.fn[X0]->res[0].begin(), t.fn[X0]->res[0].end(),
                                          x0_cache[k].begin()); }
                if (t.fn[U0]) { fill_p_only(*t.fn[U0]); t.fn[U0]->call();
                                std::copy(t.fn[U0]->res[0].begin(), t.fn[U0]->res[0].end(),
                                          u0_cache[k].begin()); }
            }
        }
    };

    // Driver. The interior-point algorithm is built ONCE (build() is where all
    // allocation happens) and reused: Fatrop's optimize() internally calls
    // reset() + initializer_->initialize(), so each call is a fresh cold solve
    // that re-reads get_initial_xk/uk and the (mutated) parameters. solve() does
    // no dynamic allocation of its own — params, initial-guess caches and the
    // x_sol/u_sol read-back buffers are all pre-sized during setup().
    struct Solver
    {
        std::shared_ptr<GeneratedOcp> ocp;
        std::map<std::string, double> opt_d;
        std::map<std::string, long long> opt_i;
        std::map<std::string, bool> opt_b;

        // Built once in setup().
        fatrop::OptionRegistry options;
        std::shared_ptr<fatrop::NlpOcp> nlp;
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
            ocp->allocate_caches();
            nlp = std::make_shared<fatrop::NlpOcp>(ocp);
            builder = std::make_unique<fatrop::IpAlgBuilder<fatrop::OcpType>>(nlp);
            ipalg = builder->with_options_registry(&options).build();
            ipdata = builder->get_ipdata();
            // Options must be applied AFTER build() (the registry is populated then).
            for (auto &kv : opt_d) options.set_option<double>(kv.first, kv.second);
            for (auto &kv : opt_i) options.set_option<fatrop::Index>(kv.first, (fatrop::Index)kv.second);
            for (auto &kv : opt_b) options.set_option<bool>(kv.first, kv.second);
            int K = ocp->get_horizon_length();
            x_sol.resize(K);
            u_sol.resize(K);
            for (int k = 0; k < K; ++k)
            {
                x_sol[k].assign(ocp->get_nx(k), 0.0);
                u_sol[k].assign(ocp->get_nu(k), 0.0);
            }
            is_setup = true;
        }

        int solve()
        {
            if (!is_setup) setup();          // one-time
            ocp->fill_initial();             // no allocation (caches pre-sized)
            ret_flag = (int)ipalg->optimize();
            iter_count = (int)ipdata->iteration_number();
            solve_time = ipdata->timing_statistics().full_algorithm.elapsed();

            const auto &x = ipalg->solution_primal();
            const auto &info = ipalg->info();
            int K = ocp->get_horizon_length();
            for (int k = 0; k < K; ++k)
            {
                int nx = ocp->get_nx(k), nu = ocp->get_nu(k);
                int ox = info.offsets_primal_x[k], ou = info.offsets_primal_u[k];
                for (int i = 0; i < nx; ++i) x_sol[k][i] = x(ox + i);
                for (int i = 0; i < nu; ++i) u_sol[k][i] = x(ou + i);
            }
            return ret_flag;
        }
    };
} // namespace ocp_direct
