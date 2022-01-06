#ifndef CASADICODEGENINCLUDED
#define CASADICODEGENINCLUDED
#include <vector>
#include <string>
#include <memory>
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
    class fatrop_eval : public fatrop_eval_base
    {
    public:
        /// constructor from file
        fatrop_eval(shared_ptr<void> handle, const std::string &function_name);
        void codegen_init(void *handle, const std::string &function_name);
        /// pointer to result_buffer
        double *output_buffer_p;
        /// casadi int work vector
        casadi_int *iw;
        /// casadi double work vector
        double *w;
        /// increase reference counter
        signal_t incref;
        /// decrease reference counter
        signal_t decref;
        /// input size
        int *input_size;
        /// release casadi memory
        casadi_release_t release;
        /// thread local mem id
        int mem;
        /// double work vector
        vector<double> work_vector_d;
        /// int work vector
        vector<casadi_int> work_vector_i;
        /// evaluate function and save res in "ccs format with lda==out_m"
        int eval(const double **arg, double *res);
        /// for reference counting of handle pointer
        shared_ptr<void> handle;
        ~fatrop_eval();
    };
} // fatrop

#endif // CASADICODEGENINCLUDED