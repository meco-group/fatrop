#ifndef BASICOCPSAMPLERSINCLUDED
#define BASICOCPSAMPLERSINCLUDED
#include "function_evaluation/FunctionEvaluation.hpp"
#include "ocp/OCPAdapter.hpp"
#include <memory>
// todo change the name of this file to evaluators
namespace fatrop
{
    class StageExpression
    {
    public:
        virtual void Eval(const double *u, const double *x, const double *global_params, const double *stage_params, double *res) = 0;
        virtual int n_rows() = 0;
        virtual int n_cols() = 0;
        int Size();
    };
    class IndexEpression : public StageExpression
    {
    public:
        IndexEpression(const bool control, const vector<int> offsets_in, const vector<int> offsets_out);
        void Eval(const double *u, const double *x, const double *global_params, const double *stage_params, double *res) override;
        int n_rows() override;
        int n_cols() override;

    private:
        const int _no_var;
        const vector<int> _offsets_in;
        const vector<int> _offsets_out;
        const bool _control;
    };
    class EvalBaseSE : public StageExpression
    {
    public:
        EvalBaseSE(const shared_ptr<EvalBase> &evalbase);
        void Eval(const double *u, const double *x, const double *global_params, const double *stage_params, double *res) override;
        int n_rows() override;
        int n_cols() override;

    private:
        shared_ptr<EvalBase> evalbase_;
        const int n_rows_;
        const int n_cols_;
    };
    class StageExpressionEvaluatorBase
    {
    public:
        virtual int Size() const = 0;
        virtual int n_rows() const= 0;
        virtual int n_cols() const = 0;
        virtual int Evaluate(const vector<double> &solution, const vector<double> &global_params, const vector<double> &stage_params, vector<double> &result) const = 0;
        vector<double> Evaluate(const vector<double> &solution, const vector<double> &global_params, const vector<double> &stage_params) const
        {
            vector<double> res(Size());
            Evaluate(solution, global_params, stage_params, res);
            return res;
        };

    protected:
        // virtual int Evaluate()
    };
    class StageControlGridSampler : public StageExpressionEvaluatorBase
    {
    public:
        StageControlGridSampler(int nu, int nx, int no_stage_params, int K, const shared_ptr<StageExpression> &eval);
        int Evaluate(const vector<double> &solution, const vector<double> &global_params, const vector<double> &stage_params, vector<double> &result) const;
        int Size() const;
        int n_rows() const;
        int n_cols() const;
        int K() const;
        const int nu;
        const int nx;
        const int no_stage_params;
        const int K_;

    private:
        shared_ptr<StageExpression> eval_;
    };

    class OCPTimeStepSampler : public StageExpressionEvaluatorBase
    {
    public:
        OCPTimeStepSampler(int nu, int nx, int no_stage_params, int K, int k, const shared_ptr<StageExpression> &eval) : nu(nu), nx(nx), no_stage_params(no_stage_params), K_(K), k_(k), eval_(eval)
        {
        }
        int Evaluate(const vector<double> &solution, const vector<double> &global_params, const vector<double> &stage_params, vector<double> &result) const
        {
            eval_->Eval(solution.data() + (k_ < K_ - 1 ? k_ * (nu + nx) : (K_ - 2) * (nu + nx)), solution.data() + (k_ < K_ - 1 ? k_ * (nu + nx) + nu : (K_ - 1) * (nu + nx)), global_params.data(), stage_params.data() + k_ * no_stage_params, result.data());
            return 0;
        }
        int Size() const
        {
            return eval_->Size();
        }
        int n_rows() const
        {
            return eval_->n_rows();
        }
        int n_cols() const
        {
            return eval_->n_cols();
        }
        const int nu;
        const int nx;
        const int no_stage_params;
        const int K_;
        const int k_;

    private:
        shared_ptr<StageExpression> eval_;
    };
    class StageExpressionEvaluatorFactory
    {
    public:
        StageExpressionEvaluatorFactory(const shared_ptr<StageExpression> &eval, int nu, int nx, int no_stage_params, int K) : nu(nu), nx(nx), no_stage_params(no_stage_params), K(K), eval_(eval){};
        // at_t0()
        OCPTimeStepSampler at_t0()
        {
            return OCPTimeStepSampler(nu, nx, no_stage_params, K, 0, eval_);
        }
        // at_tf()
        OCPTimeStepSampler at_tf()
        {
            return OCPTimeStepSampler(nu, nx, no_stage_params, K, K - 1, eval_);
        }
        // at_tk(int k)
        OCPTimeStepSampler at_tk(int k)
        {
            return OCPTimeStepSampler(nu, nx, no_stage_params, K, k, eval_);
        }
        // evaluate at control grid
        StageControlGridSampler at_control()
        {
            return StageControlGridSampler(nu, nx, no_stage_params, K, eval_);
        }
        const int nu;
        const int nx;
        const int no_stage_params;
        const int K;
        shared_ptr<StageExpression> eval_;
    };

    class ParameterSetter
    {
    public:
        ParameterSetter(const vector<int> &offsets_in, const vector<int> &offsets_out, const int no_stage_params, const int no_var, const int K, const bool global);
        void SetValue(vector<double> &global_params, vector<double> &stage_params, const double value[]);
        void SetValue(vector<double> &global_params, vector<double> &stage_params, const initializer_list<double> il_);

    private:
        const vector<int> _offsets_in;
        const vector<int> _offsets_out;
        const int no_stage_params;

    protected:
        const int _no_var;
        const int K;
        const bool _global;
    };

} // namespace fatrop

#endif //  BASICOCPSAMPLERSINCLUDED