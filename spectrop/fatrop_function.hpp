#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "casadi_utilities.hpp"
#include "ocp/StageOCPApplication.hpp"
namespace fatrop
{
    namespace spectrop
    {
        namespace cs = casadi;
        class FatropFunction : public cs::Callback
        {
        public:
            FatropFunction(const std::shared_ptr<OCPAbstract> &ocp_impl) : app(ocp_impl)
            {
                app.build();
                fatrop::OCPDims dims = app.get_ocp_dims();
                n_vars = dims.n_u_tot + dims.n_x_tot;
                n_stage_params = dims.n_stage_params_tot;
                n_global_params = dims.n_global_params;
                arg_initial_vars = std::vector<double>(n_vars);
                arg_stage_parameters = std::vector<double>(n_stage_params);
                arg_global_parameters = std::vector<double>(n_global_params);
                cs::Callback::construct("FatropFunction");
            }
            casadi_int get_n_in() override
            {
                return 3; // initial vars, stage_parameters, global_paramets
            }
            casadi_int get_n_out() override
            {
                return 1;
            }
            cs::Sparsity get_sparsity_in(casadi_int i) override
            {
                if (i == 0)
                {
                    return cs::Sparsity::dense(n_vars, 1);
                }
                else if (i == 1)
                {
                    return cs::Sparsity::dense(n_stage_params, 1);
                }
                // else
                // {
                return cs::Sparsity::dense(n_global_params, 1);
                // }
            }
            cs::Sparsity get_sparsity_out(casadi_int i) override
            {
                return cs::Sparsity::dense(n_vars, 1);
            }

            std::string get_name_in(casadi_int i) override
            {
                if (i == 0)
                {
                    return "initial_vars";
                }
                else if (i == 1)
                {
                    return "stage_parameters";
                }
                // else
                // {
                return "global_parameters";
                // }
            }

            std::string get_name_out(casadi_int i) override
            {
                return "vars";
            }

            std::vector<cs::DM> eval(const std::vector<cs::DM> &args) const override
            {
                // inefficient implementation with dynamic memory allocations
                std::vector<double> initial_vars = args[0].nonzeros();
                std::vector<double> stage_parameters = args[1].nonzeros();
                std::vector<double> global_parameters = args[2].nonzeros();
                app.set_initial(initial_vars);
                app.set_params(global_parameters, stage_parameters);
                app.optimize();
                auto retv = std::vector<double>(n_vars);
                app.last_solution_primal().copyto(retv);
                auto ret = std::vector<cs::DM>{retv};
                return ret;
            }
            int eval_buffer(const double **arg, const std::vector<casadi_int> &sizes_arg,
                            double **res, const std::vector<casadi_int> &sizes_res) const override
            {
                // no dynamic memory allocations here
                std::copy(arg[0], arg[0] + n_vars, (double *)arg_initial_vars.data());
                std::copy(arg[1], arg[1] + n_stage_params, (double *)arg_stage_parameters.data());
                std::copy(arg[2], arg[2] + n_global_params, (double *)arg_global_parameters.data());
                app.set_initial(arg_initial_vars);
                app.set_params(arg_global_parameters, arg_stage_parameters);
                app.optimize();
                double *last_sol = ((VEC *)app.last_solution_primal())->pa;
                if (res[0])
                    std::copy(last_sol, last_sol + n_vars, res[0]);
                return 0;
            }
            bool has_eval_buffer() const override { return true; };

        private:
            std::vector<double> arg_initial_vars;
            mutable std::vector<double> arg_stage_parameters;
            mutable std::vector<double> arg_global_parameters;
            fatrop::OCPApplication app;
            int n_vars;
            int n_stage_params;
            int n_global_params;
        };

    } // namespace spectrop
} // namespace fatrop