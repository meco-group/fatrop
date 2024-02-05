#pragma once
#include <vector>
#include <memory>
#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif
// c++ stuff here
// TODO: this can become part of the default fatrop API as well. It has no dependency on casadi.
// forward declaration of OCPApplication
namespace fatrop
{
  class OCPApplication;
  namespace spectool
  {

    struct C_api_userdata
    {
      // constructor
      C_api_userdata(const std::shared_ptr<fatrop::OCPApplication> &app);
      const std::shared_ptr<fatrop::OCPApplication> app_;
      std::vector<std::vector<casadi_int>> sparsity_in;
      std::vector<std::vector<casadi_int>> sparsity_out;
      std::vector<double> arg_initial_vars;
      std::vector<double> arg_stage_parameters;
      std::vector<double> arg_global_parameters;
      int n_vars;
      int n_stage_params;
      int n_global_params;
      int ref_count = 0;
    };
  }
}

#ifdef __cplusplus
extern "C"
{
#endif

/* How to prefix internal symbols */
#ifdef CASADI_CODEGEN_PREFIX
#define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)
#define _CASADI_NAMESPACE_CONCAT(NS, ID) NS##ID
#define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
#define CASADI_PREFIX(ID) f_##ID
#endif

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
#if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#if defined(STATIC_LINKED)
#define CASADI_SYMBOL_EXPORT
#else
#define CASADI_SYMBOL_EXPORT __declspec(dllexport)
#endif
#elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
#define CASADI_SYMBOL_EXPORT __attribute__((visibility("default")))
#else
#define CASADI_SYMBOL_EXPORT
#endif
#endif

  CASADI_SYMBOL_EXPORT int fatrop_func(const casadi_real **arg, casadi_real **res, casadi_int *iw, casadi_real *w, int mem, void *user_data);
  CASADI_SYMBOL_EXPORT int fatrop_func_alloc_mem(void *user_data);
  CASADI_SYMBOL_EXPORT int fatrop_func_init_mem(int mem, void *user_data);
  CASADI_SYMBOL_EXPORT void fatrop_func_free_mem(int mem, void *user_data);
  CASADI_SYMBOL_EXPORT int fatrop_func_checkout(void *user_data);
  CASADI_SYMBOL_EXPORT void fatrop_func_release(int mem, void *user_data);
  CASADI_SYMBOL_EXPORT void fatrop_func_incref(void *user_data);
  CASADI_SYMBOL_EXPORT void fatrop_func_decref(void *user_data);
  CASADI_SYMBOL_EXPORT casadi_int fatrop_func_n_in(void *user_data);
  CASADI_SYMBOL_EXPORT casadi_int fatrop_func_n_out(void *user_data);
  CASADI_SYMBOL_EXPORT casadi_real fatrop_func_default_in(casadi_int i, void *user_data);
  CASADI_SYMBOL_EXPORT const char *fatrop_func_name_in(casadi_int i, void *user_data);
  CASADI_SYMBOL_EXPORT const char *fatrop_func_name_out(casadi_int i, void *user_data);
  CASADI_SYMBOL_EXPORT const casadi_int *fatrop_func_sparsity_in(casadi_int i, void *user_data);
  CASADI_SYMBOL_EXPORT const casadi_int *fatrop_func_sparsity_out(casadi_int i, void *user_data);
  CASADI_SYMBOL_EXPORT int fatrop_func_work(casadi_int *sz_arg, casadi_int *sz_res, casadi_int *sz_iw, casadi_int *sz_w, void *user_data);

#ifdef __cplusplus
} /* extern "C" */
#endif
