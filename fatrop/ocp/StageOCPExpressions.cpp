#include "StageOCPExpressions.hpp"
using namespace fatrop;
int StageExpression::size() const
{
    return n_rows() * n_cols();
}
IndexEpression::IndexEpression(const bool control, const vector<int> offsets_in, const vector<int> offsets_out) : _no_var(offsets_in.size()),
                                                                                                                  _offsets_in(offsets_in),
                                                                                                                  _offsets_out(offsets_out),
                                                                                                                  _control(control)
{
}
void IndexEpression::evaluate(const double *u, const double *x, const double *global_params, const double *stage_params, double *res) const
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
int IndexEpression::n_rows() const
{
    return _no_var;
}
int IndexEpression::n_cols() const
{
    return 1;
}
EvalBaseSE::EvalBaseSE(const shared_ptr<EvalBase> &evalbase) : evalbase_(evalbase), n_rows_(evalbase->out_m), n_cols_(evalbase->out_n)
{
}
void EvalBaseSE::evaluate(const double *u, const double *x, const double *global_params, const double *stage_params, double *res) const
{
    const double *arg[] = {u, x, stage_params, global_params};
    evalbase_->eval_array(arg, res);
}
int EvalBaseSE::n_rows() const
{
    return n_rows_;
}
int EvalBaseSE::n_cols() const
{
    return n_cols_;
}
StageControlGridSampler::StageControlGridSampler(int nu, int nx, int no_stage_params, int K, const shared_ptr<StageExpression> &eval) : nu(nu),
                                                                                                                             nx(nx),
                                                                                                                             no_stage_params(no_stage_params),
                                                                                                                             K_(K),
                                                                                                                             eval_(eval)
{
}

int StageControlGridSampler::size() const
{
    return K_ * eval_->size();
}
int StageControlGridSampler::n_rows() const
{
    return eval_->n_rows();
}
int StageControlGridSampler::n_cols() const
{
    return eval_->n_cols();
}
int StageControlGridSampler::K() const
{
    return K_;
}

int StageControlGridSampler::evaluate(const vector<double> &solution, const vector<double> &global_params, const vector<double> &stage_params, vector<double> &sample) const
{
    const double *sol_p = solution.data();
    double *res_p = sample.data();
    const double *global_params_p = global_params.data();
    const double *stage_params_p = stage_params.data();
    int size = eval_->size();
    for (int k = 0; k < K_ - 1; k++)
    {
        eval_->evaluate(sol_p + k * (nu + nx), sol_p + k * (nu + nx) + nu, global_params_p, stage_params_p + k * no_stage_params, res_p + k * size);
    };
    eval_->evaluate(sol_p + (K_ - 2) * (nu + nx), sol_p + (K_ - 1) * (nu + nx), global_params_p, stage_params_p + (K_ - 1) * no_stage_params, res_p + (K_ - 1) * size);
    return 0;
}

ParameterSetter::ParameterSetter(const vector<int> &offsets_in, const vector<int> &offsets_out, const int no_stage_params, const int no_var, const int K, const bool global) : _offsets_in(offsets_in), _offsets_out(offsets_out), no_stage_params(no_stage_params), _no_var(no_var), K(K), _global(global)
{
}
void ParameterSetter::set_value(vector<double>& global_params, vector<double>& stage_params, const double value[])
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
void ParameterSetter::set_value(vector<double>& global_params, vector<double>& stage_params, const initializer_list<double> il_)
{
    assert((int)il_.size() == _no_var);
    set_value(global_params, stage_params, il_.begin());
}