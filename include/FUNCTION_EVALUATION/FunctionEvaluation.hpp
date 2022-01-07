#ifndef FUNCTIONEVALUATIONINCLUDED
#define FUNCTIONEVALUATIONINCLUDED
#include <vector>
#include "BLASFEO_WRAPPER/FatropLinearAlgebraBlasfeo.hpp"
using namespace std;
namespace fatrop
{
    /// Class used to evaluate a numerical functions. Functions can be implemented by hand or casadi codegen API or by plain casadi.
    class fatrop_eval_base
    {
    public:
        /// number of input vectors of the function
        int n_in;
        /// number of columns in output matrix
        int out_m;
        /// number of rows in output matrix
        int out_n;
        /// number of nonzeros in output matrix
        int out_nnz;
        /// sparsity pattern of output matrix sparsity pattern [m,n|0,ncol0, ncol0:1 , ..., | nnz | row_el0, row_el1, ...]
        vector<int> sparsity_out;
        /// buffer to safe evaluation result, in a buffer we always save a matrix in CCS format with lda==out_m
        vector<double> buffer;
        /// evaluate function and save res in "ccs format with lda==out_m"
        virtual int eval_buffer(const double **arg) = 0;
        /// evaluate function and save res in "blasfeo format"
        int eval_bf(const double **arg, MAT *bf_mat) 
        {
            double *buffer_p = buffer.data();
            // todo make this static polymorphism using CRTP
            int res = eval_buffer(arg);
            PACKMAT(out_m, out_n, buffer_p, out_m, bf_mat, 0, 0);
            return res;
        }
        ~fatrop_eval_base(){};
    };
};     // namespace fatrop
#endif // FUNCTIONEVALUATIONINCLUDED
