#include "fatrop_function_c_api.hpp"

// c++ stuff here
namespace fatrop
{
  namespace spectool
  {
    struct SparsityAux
    {
      static std::vector<casadi_int> dense(casadi_int m, casadi_int n)
      {
        std::vector<casadi_int> ret(3 + n + m * n);
        casadi_int count = 0;
        ret[count++] = m;
        ret[count++] = n;
        for (casadi_int i = 0; i < n + 1; i++)
        {
          ret[count++] = i * m;
        }
        for (casadi_int i = 0; i < n; i++)
        {
          for (casadi_int j = 0; j < m; j++)
          {
            ret[count++] = j;
          }
        }
        return ret;
      }
    };

    C_api_userdata::C_api_userdata(const std::shared_ptr<fatrop::OCPAbstract> &ocpimpl_, const cs::Dict &options) : app_(std::make_shared<fatrop::OCPApplication>(ocpimpl_))
    {
      app_->build();
      // go over the options and set
      for (auto opt : options)
      {
        if (opt.second.is_double())
          app_->set_option(opt.first, (double)opt.second);
        else if (opt.second.is_int())
          app_->set_option(opt.first, (int)opt.second);
        else if (opt.second.is_bool())
          app_->set_option(opt.first, (bool)opt.second);
      }
      fatrop::OCPDims dims = app_->get_ocp_dims();
      n_vars = dims.n_u_tot + dims.n_x_tot;
      n_stage_params = dims.n_stage_params_tot;
      n_global_params = dims.n_global_params;
      arg_initial_vars = std::vector<double>(n_vars);
      arg_stage_parameters = std::vector<double>(n_stage_params);
      arg_global_parameters = std::vector<double>(n_global_params);
      sparsity_in = {
          SparsityAux::dense(n_vars, 1),
          SparsityAux::dense(n_stage_params, 1),
          SparsityAux::dense(n_global_params, 1)};
      sparsity_out = {cs::Sparsity::dense(n_vars, 1)};
    };
  }
}

/* f:(i0)->(o0) */

CASADI_SYMBOL_EXPORT int fatrop_func(const casadi_real **arg, casadi_real **res, casadi_int *iw, casadi_real *w, int mem, void *user_data)
{
  auto *udata = static_cast<fatrop::spectool::C_api_userdata *>(user_data);
  std::copy(arg[0], arg[0] + udata->n_vars, (double *)udata->arg_initial_vars.data());
  std::copy(arg[1], arg[1] + udata->n_stage_params, (double *)udata->arg_stage_parameters.data());
  std::copy(arg[2], arg[2] + udata->n_global_params, (double *)udata->arg_global_parameters.data());
  udata->app_->set_initial(udata->arg_initial_vars);
  udata->app_->set_params(udata->arg_global_parameters, udata->arg_stage_parameters);
  int ret = udata->app_->optimize();
  if (ret != 0)
    throw std::runtime_error("fatrop solver failed");
  double *last_sol = ((VEC *)udata->app_->last_solution_primal())->pa;
  if (res[0])
    std::copy(last_sol, last_sol + udata->n_vars, res[0]);
  return 0;
}

CASADI_SYMBOL_EXPORT int fatrop_func_alloc_mem(void *user_data)
{
  return 0;
}

CASADI_SYMBOL_EXPORT int fatrop_func_init_mem(int mem, void *user_data)
{
  return 0;
}

CASADI_SYMBOL_EXPORT void fatrop_func_free_mem(int mem, void *user_data)
{
}

CASADI_SYMBOL_EXPORT int fatrop_func_checkout(void *user_data)
{
  return 0;
}

CASADI_SYMBOL_EXPORT void fatrop_func_release(int mem, void *user_data)
{
}

CASADI_SYMBOL_EXPORT void fatrop_func_incref(void *user_data)
{
  auto *udata = static_cast<fatrop::spectool::C_api_userdata *>(user_data);
  udata->ref_count++;
}

CASADI_SYMBOL_EXPORT void fatrop_func_decref(void *user_data)
{
  auto *udata = static_cast<fatrop::spectool::C_api_userdata *>(user_data);
  udata->ref_count--;
  if (udata->ref_count == 0)
  {
    delete udata;
  }
}

CASADI_SYMBOL_EXPORT casadi_int fatrop_func_n_in(void *user_data) { return 3; }

CASADI_SYMBOL_EXPORT casadi_int fatrop_func_n_out(void *user_data) { return 1; }

CASADI_SYMBOL_EXPORT casadi_real fatrop_func_default_in(casadi_int i, void *user_data)
{
  switch (i)
  {
  default:
    return 0;
  }
}

CASADI_SYMBOL_EXPORT const char *fatrop_func_name_in(casadi_int i, void *user_data)
{
  switch (i)
  {
  case 0:
    return "initial_vars";
  case 1:
    return "stage_parameters";
  case 2:
    return "global_parameters";
  default:
    return 0;
  }
}

CASADI_SYMBOL_EXPORT const char *fatrop_func_name_out(casadi_int i, void *user_data)
{
  switch (i)
  {
  case 0:
    return "solution_vars";
  default:
    return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int *fatrop_func_sparsity_in(casadi_int i, void *user_data)
{
  auto *udata = static_cast<fatrop::spectool::C_api_userdata *>(user_data);
  return udata->sparsity_in[i].data();
}

CASADI_SYMBOL_EXPORT const casadi_int *fatrop_func_sparsity_out(casadi_int i, void *user_data)
{
  auto *udata = static_cast<fatrop::spectool::C_api_userdata *>(user_data);
  return udata->sparsity_out[i].data();
}

CASADI_SYMBOL_EXPORT int fatrop_func_work(casadi_int *sz_arg, casadi_int *sz_res, casadi_int *sz_iw, casadi_int *sz_w, void *user_data)
{
  if (sz_arg)
    *sz_arg = 3;
  if (sz_res)
    *sz_res = 1;
  if (sz_iw)
    *sz_iw = 0;
  if (sz_w)
    *sz_w = 0;
  return 0;
}
