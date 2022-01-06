#ifndef CASADICODEGENINCLUDED
#define CASADICODEGENINCLUDED
#include <vector>
#include <string>
#include "FunctionEvaluation.hpp"
/* Typedefs */
typedef long long int casadi_int;
typedef void (*signal_t)(void);
typedef casadi_int (*getint_t)(void);
typedef int (*work_t)(casadi_int *sz_arg, casadi_int *sz_res, casadi_int *sz_iw, casadi_int *sz_w);
typedef const casadi_int *(*sparsity_t)(casadi_int ind);
typedef int (*eval_t)(const double **arg, double **res, casadi_int *iw, double *w, int mem);
typedef int (*casadi_checkout_t)(void);
typedef void (*casadi_release_t)(int);
namespace fatrop
{
    // class fatrop_eval : public fatrop_eval_base
    // {
    // public:
    //     /// constructor from file
    //     fatrop_eval(const fatrop_eval_construct_args<SOURCE_CASADI_CODEGEN> &constr_args, const std::string &function_name);
    //     void codegen_init(void *handle, const std::string &function_name);
    //     int *pmo_offset_p;
    //     /// pointer to result_buffer
    //     double *output_buffer_p;
    //     /// pointer to casadi codegen evalutation function
    //     eval_t eval; // !! multhithreading of this function not yet supported
    //     /// casadi int work vector
    //     casadi_int *iw;
    //     /// casadi double work vector
    //     double *w;
    //     /// increase reference counter
    //     signal_t incref;
    //     /// decrease reference counter
    //     signal_t decref;
    //     /// input size
    //     int *input_size;
    //     /// release casadi memory
    //     casadi_release_t release;
    //     /// thread local mem id
    //     int mem;
    //     /// output_buffer
    //     std::vector<double> output_buffer;
    //     /// double work vector
    //     std::vector<double> work_vector_d;
    //     /// int work vector
    //     std::vector<casadi_int> work_vector_i;

    //     /// evaluate function and save res in "panel-major-order format" this is the format used in blasfeo (panel-size 4)
    //     int eval_pmo(const double **arg, double *res);
    //     /// evaluate function and save res in "column-major-order format"
    //     int eval_cmo(const double **arg, double *res);
    //     /// evaluate function and save res in "ccs format"
    //     int eval_ccs(const double **arg, double *res);
    //     ~fatrop_eval();
    // };
} // fatrop

#endif // CASADICODEGENINCLUDED