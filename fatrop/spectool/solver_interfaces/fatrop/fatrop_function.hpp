#pragma once
#include <casadi/casadi.hpp>
#include <casadi/core/function_internal.hpp>
#include <string>
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
#include "fatrop/ocp/StageOCPApplication.hpp"
#include "fatrop/ocp/OCPDims.hpp"
#include "fatrop/ocp/CasadiCApiUserdataWrap.hpp"
namespace fatrop
{
    namespace spectool
    {
        namespace cs = casadi;
        class FatropFunctionInternal : public cs::FunctionInternal
        {
        public:
            FatropFunctionInternal(const std::string &name, const std::shared_ptr<fatrop::OCPAbstractApplication> &app_) : cs::FunctionInternal(name), app(app_)
            {
                // cs::CallbackInternal::construct(name, cb);
                fatrop::OCPDims dims = app->get_ocp_dims();
                n_vars = dims.n_u_tot + dims.n_x_tot;
                n_stage_params = dims.n_stage_params_tot;
                n_global_params = dims.n_global_params;
                arg_initial_vars = std::vector<double>(n_vars);
                arg_stage_parameters = std::vector<double>(n_stage_params);
                arg_global_parameters = std::vector<double>(n_global_params);
            }
            size_t get_n_in() override
            {
                return 3; // initial vars, stage_parameters, global_paramets
            }
            size_t get_n_out() override
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

            int eval(const double **arg,
                     double **res, casadi_int *iw, double *w, void *mem) const override
            {
                // no dynamic memory allocations here
                std::copy(arg[0], arg[0] + n_vars, (double *)arg_initial_vars.data());
                std::copy(arg[1], arg[1] + n_stage_params, (double *)arg_stage_parameters.data());
                std::copy(arg[2], arg[2] + n_global_params, (double *)arg_global_parameters.data());
                app->set_initial(arg_initial_vars);
                app->set_params(arg_global_parameters, arg_stage_parameters);
                int ret = app->optimize();
                double *last_sol = ((VEC *)app->last_solution_primal())->pa;
                if (res[0])
                    std::copy(last_sol, last_sol + n_vars, res[0]);
                if(error_on_fail_ && ret != 0)
                    throw std::runtime_error("fatrop solver failed");
                return 0;
            }
            // destructor
            ~FatropFunctionInternal() override
            {
                clear_mem();
            }

            std::string class_name() const override { return "FatropFunction"; }
            std::vector<double> arg_initial_vars;
            std::vector<double> arg_stage_parameters;
            std::vector<double> arg_global_parameters;
            std::shared_ptr<fatrop::OCPAbstractApplication> app;
            int n_vars;
            int n_stage_params;
            int n_global_params;
        };

        class FatropFunction : public cs::Function
        {
        public:
            FatropFunction(const std::string&name,  const std::shared_ptr<fatrop::OCPAbstractApplication> app_, const cs::Dict &options, const cs::Dict& funct_opts)
            {
                if (!is_null())
                {
                    casadi_error("Cannot create 'FatropFunction': Internal class already created");
                }
                auto ptr = new FatropFunctionInternal(std::string("FatropFunction_") + name, app_);
                static_cast<cs::SharedObject *>(this)->own(ptr);
                auto options_ = cs::Dict();
                ptr -> error_on_fail_ = cs::get_from_dict(funct_opts, "error_on_fail", true);
                static_cast<cs::Function *>(this)->operator->()->casadi::FunctionInternal::construct(options_);
            }

        private:
        };

    } // namespace spectool
} // namespace fatrop