#pragma once
#include <casadi/casadi.hpp>
#include <casadi/core/callback_internal.hpp>
#include <string>
#include "casadi_utilities.hpp"
#include "ocp/StageOCPApplication.hpp"
namespace fatrop
{
    namespace spectrop
    {
        namespace cs = casadi;
        class FatropFunctionInternal : public cs::CallbackInternal
        {
        public:
            FatropFunctionInternal(const std::string &name, cs::Callback *cb, const std::shared_ptr<fatrop::OCPApplication> &app_) : cs::CallbackInternal(name, cb), app(app_)
            {
                // cs::CallbackInternal::construct(name, cb);
            }
            std::vector<double> arg_initial_vars;
            std::vector<double> arg_stage_parameters;
            std::vector<double> arg_global_parameters;
            std::shared_ptr<fatrop::OCPApplication> app;
            int n_vars;
            int n_stage_params;
            int n_global_params;
        };
        class FatropFunction : public cs::Callback
        {
        public:
            FatropFunction(const std::shared_ptr<fatrop::OCPAbstract> &ocpimpl_)
            {
                auto app_ = std::make_shared<fatrop::OCPApplication>(ocpimpl_);
                app_->build();
                // cs::Callback::construct("FatropFunction");
                if (!is_null())
                {
                    casadi_error("Cannot create 'FatropFunction': Internal class already created");
                }
                auto ptr = new FatropFunctionInternal("FatropFunction", this, app_);
                static_cast<cs::SharedObject *>(this)->own(ptr);
                auto options_ = cs::Dict();
                auto internal_p = static_cast<FatropFunctionInternal *>(this->get());
                fatrop::OCPDims dims = internal_p->app->get_ocp_dims();
                internal_p->n_vars = dims.n_u_tot + dims.n_x_tot;
                internal_p->n_stage_params = dims.n_stage_params_tot;
                internal_p->n_global_params = dims.n_global_params;
                internal_p->arg_initial_vars = std::vector<double>(internal_p->n_vars);
                internal_p->arg_stage_parameters = std::vector<double>(internal_p->n_stage_params);
                internal_p->arg_global_parameters = std::vector<double>(internal_p->n_global_params);
                static_cast<cs::Function *>(this)->operator->()->casadi::FunctionInternal::construct(options_);
            }
            ~FatropFunction() override
            {
                // bypass the destructor of cs::Callback
                // this -> cs::Function::~Function();
            };
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
                auto internal_p = static_cast<FatropFunctionInternal *>(this->get());
                if (i == 0)
                {
                    return cs::Sparsity::dense(internal_p->n_vars, 1);
                }
                else if (i == 1)
                {
                    return cs::Sparsity::dense(internal_p->n_stage_params, 1);
                }
                // else
                // {
                return cs::Sparsity::dense(internal_p->n_global_params, 1);
                // }
            }
            cs::Sparsity get_sparsity_out(casadi_int i) override
            {
                auto internal_p = static_cast<FatropFunctionInternal *>(this->get());
                return cs::Sparsity::dense(internal_p->n_vars, 1);
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
                auto internal_p = static_cast<FatropFunctionInternal *>(this->get());
                // inefficient implementation with dynamic memory allocations
                std::vector<double> initial_vars = args[0].nonzeros();
                std::vector<double> stage_parameters = args[1].nonzeros();
                std::vector<double> global_parameters = args[2].nonzeros();
                internal_p->app->set_initial(initial_vars);
                internal_p->app->set_params(global_parameters, stage_parameters);
                internal_p->app->optimize();
                auto retv = std::vector<double>(internal_p->n_vars);
                internal_p->app->last_solution_primal().copyto(retv);
                auto ret = std::vector<cs::DM>{retv};
                return ret;
            }
            int eval_buffer(const double **arg, const std::vector<casadi_int> &sizes_arg,
                            double **res, const std::vector<casadi_int> &sizes_res) const override
            {
                // no dynamic memory allocations here
                auto internal_p = static_cast<FatropFunctionInternal *>(this->get());
                std::copy(arg[0], arg[0] + internal_p->n_vars, (double *)internal_p->arg_initial_vars.data());
                std::copy(arg[1], arg[1] + internal_p->n_stage_params, (double *)internal_p->arg_stage_parameters.data());
                std::copy(arg[2], arg[2] + internal_p->n_global_params, (double *)internal_p->arg_global_parameters.data());
                internal_p->app->set_initial(internal_p->arg_initial_vars);
                internal_p->app->set_params(internal_p->arg_global_parameters, internal_p->arg_stage_parameters);
                internal_p->app->optimize();
                double *last_sol = ((VEC *)internal_p->app->last_solution_primal())->pa;
                if (res[0])
                    std::copy(last_sol, last_sol + internal_p->n_vars, res[0]);
                return 0;
            }
            bool has_eval_buffer() const override { return true; };

        private:
        };

    } // namespace spectrop
} // namespace fatrop