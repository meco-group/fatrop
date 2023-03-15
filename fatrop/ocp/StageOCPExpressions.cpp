#include "StageOCPExpressions.hpp"
using namespace fatrop;
int StageExpression::Size()
{
    return n_rows() * n_cols();
}
IndexEpression::IndexEpression(const bool control, const vector<int> offsets_in, const vector<int> offsets_out) : _no_var(offsets_in.size()),
                                                                                                                  _offsets_in(offsets_in),
                                                                                                                  _offsets_out(offsets_out),
                                                                                                                  _control(control)
{
}
void IndexEpression::Eval(const double *u, const double *x, const double *global_params, const double *stage_params, double *res)
{
    if (_control)
    {
        for (int i = 0; i < _no_var; i++)
        {
            res[_offsets_in.at(i)] = u[_offsets_out.at(i)];
        }
    }
    else
    {
        for (int i = 0; i < _no_var; i++)
        {
            res[_offsets_in.at(i)] = x[_offsets_out.at(i)];
        }
    }
};
int IndexEpression::n_rows()
{
    return _no_var;
}
int IndexEpression::n_cols()
{
    return 1;
}
EvalBaseSE::EvalBaseSE(const shared_ptr<EvalBase> &evalbase) : evalbase_(evalbase), n_rows_(evalbase->out_m), n_cols_(evalbase->out_n)
{
}
void EvalBaseSE::Eval(const double *u, const double *x, const double *global_params, const double *stage_params, double *res)
{
    const double *arg[] = {u, x, stage_params, global_params};
    evalbase_->eval_array(arg, res);
}
int EvalBaseSE::n_rows()
{
    return n_rows_;
}
int EvalBaseSE::n_cols()
{
    return n_cols_;
}
StageOCPControlSampler::StageOCPControlSampler(int nu, int nx, int no_stage_params, int K, const shared_ptr<StageExpression> &eval) : nu(nu),
                                                                                                                             nx(nx),
                                                                                                                             no_stage_params(no_stage_params),
                                                                                                                             K_(K),
                                                                                                                             eval_(eval)
{
}

int StageOCPControlSampler::Size()
{
    return K_ * eval_->Size();
}
int StageOCPControlSampler::n_rows()
{
    return eval_->n_rows();
}
int StageOCPControlSampler::n_cols()
{
    return eval_->n_cols();
}
int StageOCPControlSampler::K()
{
    return K_;
}

int StageOCPControlSampler::Evaluate(const vector<double> &solution, const vector<double> &global_params, const vector<double> &stage_params, vector<double> &sample)
{
    const double *sol_p = solution.data();
    double *res_p = sample.data();
    const double *global_params_p = global_params.data();
    const double *stage_params_p = stage_params.data();
    int size = eval_->Size();
    for (int k = 0; k < K_ - 1; k++)
    {
        eval_->Eval(sol_p + k * (nu + nx), sol_p + k * (nu + nx) + nu, global_params_p, stage_params_p + k * no_stage_params, res_p + k * size);
    };
    eval_->Eval(sol_p + (K_ - 2) * (nu + nx), sol_p + (K_ - 1) * (nu + nx), global_params_p, stage_params_p + (K_ - 1) * no_stage_params, res_p + (K_ - 1) * size);
    return 0;
}

ParameterSetter::ParameterSetter(const vector<int> &offsets_in, const vector<int> &offsets_out, const int no_stage_params, const int no_var, const int K, const bool global) : _offsets_in(offsets_in), _offsets_out(offsets_out), no_stage_params(no_stage_params), _no_var(no_var), K(K), _global(global)
{
}
void ParameterSetter::SetValue(vector<double>& global_params, vector<double>& stage_params, const double value[])
{
    if (_global)
    {
        double *params = global_params.data();
        for (int i = 0; i < _no_var; i++)
        {
            params[_offsets_out.at(i)] = value[_offsets_in.at(i)];
        }
    }
    else // stage paramter
    {
        double *params = stage_params.data(); 
        for (int k = 0; k < K; k++)
        {
            for (int i = 0; i < _no_var; i++)
            {
                params[_offsets_out.at(i) + k * no_stage_params] = value[_offsets_in.at(i)];
            }
        }
    }
}
void ParameterSetter::SetValue(vector<double>& global_params, vector<double>& stage_params, const initializer_list<double> il_)
{
    assert((int)il_.size() == _no_var);
    SetValue(global_params, stage_params, il_.begin());
}