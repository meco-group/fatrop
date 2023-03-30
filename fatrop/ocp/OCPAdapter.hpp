#ifndef OCPEVALUATORINCLUDED
#define OCPEVALUATORINCLUDED
#include "OCPKKT.hpp"
#include "OCPAbstract.hpp"
#include "solver/FatropData.hpp"
#include "OCP.hpp"
#include <memory>
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP->name)
#define AUXMACRO(type, name, suffix) type name##suffix = ((type)OCP->aux.name)
#define SOLVERMACRO(type, name, suffix) type name##suffix = ((type)name)
// #ifdef ENABLE_MULTITHREADING
// #include <omp.h>
// #endif

// Example: you can make a for loop parallel with the following code. But to be further investigated to do this in an efficient way...
// #ifdef ENABLE_MULTITHREADING
// #pragma omp parallel for
// <your for loop>
// #endif

namespace fatrop
{
    class OCPAdapter : public OCP // public OCP -> also include KKTmemory, OCPDims, ...
    {
    public:
        OCPAdapter(const std::shared_ptr<OCPAbstract> &ocptempl_) : K(ocptempl_->get_horizon_length()),
                                                           nuexpr(TransformRange<int>(0, K, [&ocptempl_](int k)
                                                                                      { return ocptempl_->get_nuk(k); })),
                                                           nxexpr(TransformRange<int>(0, K, [&ocptempl_](int k)
                                                                                      { return ocptempl_->get_nxk(k); })),
                                                           ngexpr(TransformRange<int>(0, K, [&ocptempl_](int k)
                                                                                      { return ocptempl_->get_ngk(k); })),
                                                           ngineqexpr(TransformRange<int>(0, K, [&ocptempl_](int k)
                                                                                          { return ocptempl_->get_ng_ineq_k(k); })),
                                                           nstageparamsexpr(TransformRange<int>(0, K, [&ocptempl_](int k)
                                                                                                { return ocptempl_->get_n_stage_params_k(k); })),
                                                           offs_stageparams(offsets(nstageparamsexpr)), stageparams(sum(nstageparamsexpr), 0.0), globalparams(ocptempl_->get_n_global_params(), 0.0), ocptempl(ocptempl_)
        {
            // initialize the default parameters
            ocptempl_->get_default_global_params(globalparams.data());
            int offs = 0;
            for (int k = 0; k < K; k++)
            {
                ocptempl_->get_default_stage_paramsk(stageparams.data() + offs, k);
                offs += ocptempl_->get_n_stage_params_k(k);
            }
            x_dummy = std::vector<double>(max(nxexpr), 0.0);
        }
        int eval_lag_hess(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            const FatropVecBF &lam) override;
        int eval_constr_jac(
            OCPKKTMemory *OCP,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars) override;
        int eval_contr_viol(
            OCPKKTMemory *OCP,
            const FatropVecBF &primal_vars,
            const FatropVecBF &slack_vars,
            FatropVecBF &constraint_violation) override;
        int eval_obj_grad(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            FatropVecBF &gradient) override;
        int eval_obj(
            OCPKKTMemory *OCP,
            double obj_scale,
            const FatropVecBF &primal_vars,
            double &res);
        int integrate_dynamics(
            OCPKKTMemory *OCP,
            const int k,
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
        std::vector<double> & get_global_parameters_vec()
        {
            return globalparams;
        }
        std::vector<double> & get_stage_parameters_vec()
        {
            return stageparams;
        }
        int get_bounds(
            FatropVecBF &lower,
            FatropVecBF &upper) const override
        {
            int offs = 0;
            double *lower_p = ((VEC *)lower)->pa;
            double *upper_p = ((VEC *)upper)->pa;
            for (int k = 0; k < K; k++)
            {
                ocptempl->get_boundsk(lower_p + offs, upper_p + offs, k);
                offs += ocptempl->get_ng_ineq_k(k);
            }
            return 0;
        };
        int get_initial_sol_guess(
            FatropVecBF &initial) const override
        {
            int offs = 0;
            for (int k = 0; k < K-1; k++)
            {
                ocptempl->get_initial_uk(((VEC *)initial)->pa + offs, k);
                offs += ocptempl->get_nuk(k);
                ocptempl->get_initial_xk(((VEC *)initial)->pa + offs, k);
                offs += ocptempl->get_nxk(k);
            }
            ocptempl->get_initial_xk(((VEC *)initial)->pa + offs, K-1);
            return 0;
        }
        // virtual int GetDefaultParams(
        //     FatropOptions &params)
        // {
        // };

    public:
        int K;
        FatropVector<int> nuexpr;
        FatropVector<int> nxexpr;
        FatropVector<int> ngexpr;
        FatropVector<int> ngineqexpr;
        FatropVector<int> nstageparamsexpr;
        FatropVector<int> offs_stageparams;
        std::vector<double> stageparams;
        std::vector<double> globalparams;
        std::vector<double> x_dummy;

    private:
        std::shared_ptr<OCPAbstract> ocptempl;
    };
} // namespace fatrop

#endif // OCPEVALUATORINCLUDED