#ifndef BASICOCPAPPLICATIONINCLUDED
#define BASICOCPAPPLICATIONINCLUDED
#include "ocp/StageOCP.hpp"
#include "solver/AlgBuilder.hpp"
#include "ocp/OCPAdapter.hpp"
#include "solver/FatropAlg.hpp"
#include "ocp/FatropOCP.hpp"
#include "ocp/FatropOCPBuilder.hpp"
#include "ocp/OCPAbstract.hpp"
#include "ocp/StageOCPExpressions.hpp"
#include "solver/FatropStats.hpp"
#include "ocp/OCPSolution.hpp"
#include <map>
#include "json/json.h"
#include <fstream>
#include <sstream>
#include <string>

namespace fatrop
{
    class FatropSolution
    {
    public:
        // void GetPrimalSolution(vector<double> &result);
        // defautl copy constructor
        FatropSolution(const FatropSolution &other) = default;
        const vector<double> &primal_solution() { return sol_primal_; };

    protected:
        FatropSolution();
        void set_dims(const NLPDims &dims);
        void set_primal_solution(const FatropVecBF &sol);
        void set_solution(const FatropVecBF &sol_primal, const FatropVecBF &sol_dual, const FatropVecBF &sol_zL, const FatropVecBF &sol_zU);

    protected:
        vector<double> sol_primal_;
        vector<double> sol_dual_;
        vector<double> sol_zL_;
        vector<double> sol_zU_;
        friend class NLPApplication;
    };
    // TODO move this class to a separate file
    class NLPApplication
    {
    public:
        NLPApplication();

    protected:
        void build(const shared_ptr<FatropNLP> &nlp);

    public:
        int optimize();
        // TODO: make this protected and use last_solution instead and choose other name
        const FatropVecBF &last_solution_primal() const;
        const FatropVecBF &last_solution_dual() const;
        const FatropVecBF &last_solution_zL() const;
        const FatropVecBF &last_solution_zU() const;
        FatropVecBF &initial_guess_primal() const;
        FatropVecBF &initial_guess_dual() const;
        FatropVecBF &initial_guess_zL() const;
        FatropVecBF &initial_guess_zU() const;
        FatropStats get_stats() const;
        NLPDims get_nlp_dims();
        template <typename T>
        void set_option(const string &option_name, T value);

    public:
        void set_initial(const FatropSolution &initial_guess);
        const FatropOptions &get_options() const;

    protected:
        shared_ptr<FatropOptions> fatropoptions_;
        shared_ptr<FatropData> fatropdata_;
        shared_ptr<FatropNLP> nlp_;
        bool dirty = true;

    private:
        const shared_ptr<OCPAbstract> ocp_;
        shared_ptr<Journaller> journaller_;
        shared_ptr<FatropAlg> fatropalg_;
    };

    // TODO move this class to a separate file
    class OCPApplication : public NLPApplication
    {
    public:
        OCPApplication(const shared_ptr<OCPAbstract> &ocp);

    protected:
        void build();

    public:
        using NLPApplication::set_initial;
        void set_initial(vector<double> &initial_u, vector<double> &initial_x);
        OCPDims get_ocp_dims();
        void set_params(const vector<double> &global_params, const vector<double> &stage_params)
        {
            adapter->SetParams(global_params, stage_params);
        }

    private:
        const shared_ptr<OCPAbstract> ocp_;

    protected:
        vector<double> &global_parameters();
        vector<double> &stage_parameters();
        shared_ptr<OCPAdapter> adapter;
    };

    struct StageOCPSolution : public FatropSolution
    {
    public:
        StageOCPSolution(const shared_ptr<OCP> &app);

    protected:
        StageOCPSolution();
        void set_dims(const OCPDims &dims);
        void set_parameters(const vector<double> &global_params, const vector<double> &stage_params);
        int nx;
        int nu;
        int n_stage_params;
        int n_global_params;
        int K;

    public:
        // todo make this deprecated, only use Eval
        void sample(const StageControlGridSampler &sampler, vector<double> &result) const;
        vector<double> evaluate(const StageExpressionEvaluatorBase &evaluator) const;
        void evaluate(const StageExpressionEvaluatorBase &evaluator, vector<double> &result) const;

    protected:
        vector<double> global_params;
        vector<double> stage_params;
        friend class StageOCPApplication;
    };

    // struct StageOCPApplicationAbstract : public OCPApplication
    // {
    //     StageOCPApplicationAbstract(const shared_ptr<StageOCP> &ocp);
    // };
    class StageOCPApplication : public OCPApplication
    {
    public:
        StageOCPApplication(const shared_ptr<StageOCP> &ocp);

    public:
        class AppParameterSetter : public ParameterSetter
        {
        public:
            AppParameterSetter(const shared_ptr<OCPAdapter> &adapter, const shared_ptr<ParameterSetter> &ps);

        public:
            void set_value(const double value[]);
            void set_value(const initializer_list<double> il_);

        private:
            const shared_ptr<OCPAdapter> adapter_;
        };

    public:
        StageExpressionEvaluatorFactory get_expression(const string &sampler_name);
        StageExpressionEvaluatorFactory get_expression_evaluator(const shared_ptr<StageExpression> &expr);
        AppParameterSetter get_parameter_setter(const string &setter_name);
        void build();
        int optimize();
        const StageOCPSolution &last_stageocp_solution() const;
        const FatropSolution &last_solution() const { return last_stageocp_solution(); };

    public:
        const int nx_;
        const int nu_;
        const int n_stage_params_;
        const int K_;

    private:
        StageOCPSolution last_solution_;

    protected:
        map<string, shared_ptr<StageExpression>> stage_expressions;
        map<string, shared_ptr<ParameterSetter>> param_setters;
        friend class StageOCPApplicationBuilder;
    };
    class StageOCPApplicationBuilder
    {
    public:
        static StageOCPApplication from_rockit_interface(const string &functions, const string &json_spec_file);
    };
}
#endif // BASICOCPAPPLICATIONINCLUDED