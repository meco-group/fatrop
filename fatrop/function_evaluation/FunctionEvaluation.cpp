#include "function_evaluation/FunctionEvaluation.hpp"
using namespace fatrop;
int EvalBase::eval_bf(const double **arg, MAT *bf_mat)
{
#if DEBUG
    assert(bf_mat->m >= out_m);
    assert(bf_mat->n >= out_n);
#endif
    double *buffer_p = buffer.data();
    // todo make this static polymorphism using CRTP
    int res = eval_buffer(arg);
    PACKMAT(out_m, out_n, buffer_p, out_m, bf_mat, 0, 0);
    return res;
}
int EvalBase::eval_array(const double **arg, double *array)
{
    double *buffer_p = buffer.data();
    // todo make this static polymorphism using CRTP
    int res = eval_buffer(arg);
    memcpy(array, buffer_p, out_nnz * sizeof(double));
    return res;
}