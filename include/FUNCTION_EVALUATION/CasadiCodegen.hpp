#ifndef CASADICODEGENINCLUDED
#define CASADICODEGENINCLUDED
#include <vector>
#include <string>
#include <memory>
#include "FunctionEvaluation.hpp"
#include "AUX/DynamicLib.hpp"
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
    class fatrop_eval_CasGen : public fatrop_eval_base
    {
    public:
        /// constructor from file
        fatrop_eval_CasGen(const shared_ptr<DL_loader> &handle, const std::string &function_name);
        /// pointer to result_buffer
        double *output_buffer_p;
        /// pointer to casadi codegen evalutation function
        eval_t eval; // !! multhithreading of this function not yet supported
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
        int eval_buffer(const double **arg, double *res);
        /// for reference counting of handle pointer
        shared_ptr<DL_loader> handle;
        ~fatrop_eval_CasGen();
    };
} // fatrop

#endif // CASADICODEGENINCLUDED