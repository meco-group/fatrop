#ifndef BASICOCPAPPLICATIONINCLUDED
#define BASICOCPAPPLICATIONINCLUDED
#include "ocp/BFOCPBasic.hpp"
#include "solver/AlgBuilder.hpp"
#include "ocp/BFOCPAdapter.hpp"
#include "solver/FatropAlg.hpp"
#include "ocp/FatropOCP.hpp"
#include "ocp/FatropOCPBuilder.hpp"
#include "ocp/OCPAbstract.hpp"
#include "ocp/BasicOCPSamplers.hpp"
#include "solver/FatropStats.hpp"
#include "ocp/OCPSolution.hpp"
#include <map>
#include "json/json.h"
#include <fstream>
#include <sstream>
#include <string>

namespace fatrop
{
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
        FatropVecBF &LastSolution();
        FatropVecBF &InitialGuessPrimal();
        FatropStats GetStats();
        NLPDims GetNLPDims();

    protected:
        shared_ptr<FatropParams> fatropparams_;
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
        void SetInitial(vector<double> &initial_u, vector<double> &initial_x);
        OCPDims GetOCPDims();

    private:
        const shared_ptr<OCPAbstract> ocp_;

    protected:
        shared_ptr<BFOCPAdapter> adapter;
    };
    class FatropSolution
    {
    public:
        void GetPrimalSolution(vector<double> &result);
        // defautl copy constructor
        FatropSolution(const FatropSolution &other) = default;

    protected:
        FatropSolution();
        void SetDims(const NLPDims &dims);
        void SetPrimalSolution(const FatropVecBF &sol);

    protected:
        vector<double> sol_primal;
    };
    void FatropSolution::GetPrimalSolution(vector<double> &result)
    {
        result = sol_primal;
    }

    struct BasicOCPApplicationAbstract : public OCPApplication
    {
        BasicOCPApplicationAbstract(const shared_ptr<BasicOCP> &ocp);
    };

    struct BasicOCPSolution : public FatropSolution
    {
    public:
        BasicOCPSolution(const shared_ptr<BasicOCPApplicationAbstract> &app);

    protected:
        BasicOCPSolution();
        void SetDims(const OCPDims &dims);
        void Set(const FatropVecBF &sol, const vector<double> &global_params, const vector<double> &stage_params);
        int nx;
        int nu;
        int n_stage_params;
        int n_global_params;
        int K;

    public:
        // todo make this deprecated, only use Eval
        void Sample(const shared_ptr<OCPControlSampler> &sampler, vector<double> &result);
        vector<double> Eval(const shared_ptr<BasicOCPEvaluatorBase> &evaluator) const;
        void Eval(const shared_ptr<BasicOCPEvaluatorBase> &evaluator, vector<double> &result) const;

    protected:
        vector<double> global_params;
        vector<double> stage_params;
        friend class BasicOCPApplication;
    };

    class BasicOCPApplication : public BasicOCPApplicationAbstract
    {
    public:
        BasicOCPApplication(const shared_ptr<BasicOCP> &ocp);

    public:
        class AppParameterSetter : public ParameterSetter
        {
        public:
            AppParameterSetter(const shared_ptr<BFOCPAdapter> &adapter, const shared_ptr<ParameterSetter> &ps);

        public:
            void SetValue(const double value[]);

        private:
            const shared_ptr<BFOCPAdapter> adapter_;
        };

    public:
        shared_ptr<BasicOCPEvaluatorFactory> GetEvaluator(const string &sampler_name);
        shared_ptr<AppParameterSetter> GetParameterSetter(const string &setter_name);
        void Build();
        int Optimize();
        const BasicOCPSolution &LastBasicOCPSolution();

    public:
        const int nx_;
        const int nu_;
        const int K_;

    private:
        BasicOCPSolution last_solution;

    protected:
        map<string, shared_ptr<BasicOCPEvaluatorFactory>> eval_factories;
        map<string, shared_ptr<ParameterSetter>> param_setters;
        friend class BasicOCPApplicationBuilder;
    };
    class BasicOCPApplicationBuilder
    {
    public:
        static shared_ptr<BasicOCPApplication> FromRockitInterface(const string &functions, const string &json_spec_file);
    };
}
#endif // BASICOCPAPPLICATIONINCLUDED