#ifndef BASICOCPSAMPLERSINCLUDED
#define BASICOCPSAMPLERSINCLUDED
#include "function_evaluation/FunctionEvaluation.hpp"
#include "ocp/BFOCPAdapter.hpp"
#include <memory>
// todo change the name of this file to evaluators
namespace fatrop
{
    class StageEvaluator
    {
    public:
        virtual void Eval(const double *u, const double *x, const double *global_params, const double *stage_params, double *res) = 0;
        virtual int n_rows() = 0;
        virtual int n_cols() = 0;
        int Size();
    };
    class IndexEvaluator : public StageEvaluator
    {
    public:
        IndexEvaluator(const bool control, const vector<int> offsets_in, const vector<int> offsets_out);
        void Eval(const double *u, const double *x, const double *global_params, const double *stage_params, double *res) override;
        int n_rows() override;
        int n_cols() override;

    private:
        const int _no_var;
        const vector<int> _offsets_in;
        const vector<int> _offsets_out;
        const bool _control;
    };
    class EvalBaseSE : public StageEvaluator
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
    class BasicOCPEvaluatorBase
    {
    public:
        virtual int Size() = 0;
        virtual int n_rows() = 0;
        virtual int n_cols() = 0;
        virtual int Evaluate(const vector<double> &solution, const vector<double> &global_params, const vector<double> &stage_params, vector<double> &result) = 0;
        vector<double> Evaluate(const vector<double> &solution, const vector<double> &global_params, const vector<double> &stage_params)
        {
            vector<double> res(Size());
            Evaluate(solution, global_params, stage_params, res);
            return res;
        };

    protected:
        // virtual int Evaluate()
    };
    class OCPControlSampler : public BasicOCPEvaluatorBase
    {
    public:
        OCPControlSampler(int nu, int nx, int no_stage_params, int K, const shared_ptr<StageEvaluator> &eval);
        int Evaluate(const vector<double> &solution, const vector<double> &global_params, const vector<double> &stage_params, vector<double> &result);
        int Size();
        int n_rows();
        int n_cols();
        int K();
        const int nu;
        const int nx;
        const int no_stage_params;
        const int K_;

    private:
        shared_ptr<StageEvaluator> eval_;
    };
    class BasicOCPEvaluatorFactory
    {
    public:
        BasicOCPEvaluatorFactory(const shared_ptr<StageEvaluator> &eval, int nu, int nx, int no_stage_params, int K) : nu(nu), nx(nx), no_stage_params(no_stage_params), K(K), eval_(eval){};
        // at_t0()
        // at_tf()
        shared_ptr<OCPControlSampler> at_control()
        {
            return make_shared<OCPControlSampler>(nu, nx, no_stage_params, K, eval_);
        }
        const int nu;
        const int nx;
        const int no_stage_params;
        const int K;
        shared_ptr<StageEvaluator> eval_;
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
        const int _no_var;
        const int K;
        const bool _global;
    };

} // namespace fatrop

#endif //  BASICOCPSAMPLERSINCLUDED