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
        const vector<double> &PrimalSolution() { return sol_primal_; };

    protected:
        FatropSolution();
        void SetDims(const NLPDims &dims);
        void SetPrimalSolution(const FatropVecBF &sol);
        void SetSolution(const FatropVecBF &sol_primal, const FatropVecBF &sol_dual, const FatropVecBF &sol_zL, const FatropVecBF &sol_zU);

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
        void Build(const shared_ptr<FatropNLP> &nlp);

    public:
        int Optimize();
        // TODO: make this protected and use last_solution instead and choose other name
        FatropVecBF &LastSolutionPrimal();
        FatropVecBF &LastSolutionDual();
        FatropVecBF &LastSolutionZL();
        FatropVecBF &LastSolutionZU();
        FatropVecBF &InitialGuessPrimal();
        FatropVecBF &InitialGuessDual();
        FatropVecBF &InitialGuessZL();
        FatropVecBF &InitialGuessZU();
        FatropStats GetStats();
        NLPDims GetNLPDims();
        template <typename T>
        void SetOption(const string &option_name, T value);

    public:
        void SetInitial(const FatropSolution &initial_guess);
        const FatropOptions &GetOptions() const;

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
        void Build();

    public:
        vector<double> &GlobalParameters();
        vector<double> &StageParameters();
        using NLPApplication::SetInitial;
        void SetInitial(vector<double> &initial_u, vector<double> &initial_x);
        OCPDims GetOCPDims();

    private:
        const shared_ptr<OCPAbstract> ocp_;

    protected:
        shared_ptr<OCPAdapter> adapter;
    };


    struct StageOCPSolution : public FatropSolution
    {
    public:
        StageOCPSolution(const shared_ptr<OCP> &app);

    protected:
        StageOCPSolution();
        void SetDims(const OCPDims &dims);
        void SetParams(const vector<double> &global_params, const vector<double> &stage_params);
        int nx;
        int nu;
        int n_stage_params;
        int n_global_params;
        int K;

    public:
        // todo make this deprecated, only use Eval
        void Sample(const shared_ptr<StageControlGridSampler> &sampler, vector<double> &result);
        vector<double> Eval(const shared_ptr<StageExpressionEvaluatorBase> &evaluator) const;
        void Eval(const shared_ptr<StageExpressionEvaluatorBase> &evaluator, vector<double> &result) const;

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
            void SetValue(const double value[]);
            void SetValue(const initializer_list<double> il_);

        private:
            const shared_ptr<OCPAdapter> adapter_;
        };

    public:
        shared_ptr<StageExpressionEvaluatorFactory> GetExpression(const string &sampler_name);
        shared_ptr<StageExpressionEvaluatorFactory> GetExprEvaluator(const shared_ptr<StageExpression> &expr);
        shared_ptr<AppParameterSetter> GetParameterSetter(const string &setter_name);
        void Build();
        int Optimize();
        const StageOCPSolution &LastStageOCPSolution();
        const FatropSolution &LastSolution(){return LastStageOCPSolution();};

    public:
        const int nx_;
        const int nu_;
        const int n_stage_params_;
        const int K_;

    private:
        StageOCPSolution last_solution;

    protected:
        map<string, shared_ptr<StageExpression>> stage_expressions;
        map<string, shared_ptr<ParameterSetter>> param_setters;
        friend class StageOCPApplicationBuilder;
    };
    class StageOCPApplicationBuilder
    {
    public:
        static shared_ptr<StageOCPApplication> FromRockitInterface(const string &functions, const string &json_spec_file);
    };
}
#endif // BASICOCPAPPLICATIONINCLUDED