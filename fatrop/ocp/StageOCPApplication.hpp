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
        const vector<double>& PrimalSolution(){return sol_primal_;};

    protected:
        FatropSolution();
        void SetDims(const NLPDims &dims);
        void SetPrimalSolution(const FatropVecBF &sol);
        void SetSolution(const FatropVecBF &sol_primal, const FatropVecBF &sol_dual, const FatropVecBF &sol_zL, const FatropVecBF &sol_zU)
        {
            sol_primal.copyto(sol_primal_);
            sol_dual.copyto(sol_dual_);
            sol_zL.copyto(sol_zL_);
            sol_zU.copyto(sol_zU_);
        };

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
        FatropVecBF &LastSolutionDual()
        {
            return fatropdata_->lam_curr;
        }
        FatropVecBF &LastSolutionZL()
        {
            return fatropdata_->zL_curr;
        }
        FatropVecBF &LastSolutionZU()
        {
            return fatropdata_->zU_curr;
        }
        FatropVecBF &InitialGuessPrimal();
        FatropVecBF &InitialGuessDual()
        {
            return fatropdata_->lam_init;
        }
        FatropVecBF &InitialGuessZL()
        {
            return fatropdata_->zL_init;
        }
        FatropVecBF &InitialGuessZU()
        {
            return fatropdata_->zU_init;
        }
        FatropStats GetStats();
        NLPDims GetNLPDims();
        template <typename T>
        void SetOption(const string &option_name, T value)
        {
            fatropparams_->SetOption(option_name, value);
        }

    public:
        void SetInitial(const FatropSolution &initial_guess)
        {
            InitialGuessPrimal() = initial_guess.sol_primal_;
            InitialGuessDual() = initial_guess.sol_dual_;
            InitialGuessZL() = initial_guess.sol_zL_;
            InitialGuessZU() = initial_guess.sol_zU_;
        }

    protected:
        shared_ptr<FatropOptions> fatropparams_;
        shared_ptr<FatropData> fatropdata_;
        shared_ptr<FatropNLP> nlp_;
        bool dirty = true;

    private:
        const shared_ptr<OCPAbstract> ocp_;
        shared_ptr<Journaller> journaller_;
        shared_ptr<FatropAlg> fatropalg_;
    };
    template void NLPApplication::SetOption(const string &option_name, double value);
    template void NLPApplication::SetOption(const string &option_name, int value);
    template void NLPApplication::SetOption(const string &option_name, bool value);

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
    // void FatropSolution::GetPrimalSolution(vector<double> &result)
    // {
    //     result = sol_primal_;
    // }

    struct StageOCPApplicationAbstract : public OCPApplication
    {
        StageOCPApplicationAbstract(const shared_ptr<StageOCP> &ocp);
    };

    struct StageOCPSolution : public FatropSolution
    {
    public:
        StageOCPSolution(const shared_ptr<StageOCPApplicationAbstract> &app);

    protected:
        StageOCPSolution();
        void SetDims(const OCPDims &dims);
        void Set(const FatropVecBF &sol, const vector<double> &global_params, const vector<double> &stage_params);
        int nx;
        int nu;
        int n_stage_params;
        int n_global_params;
        int K;

    public:
        // todo make this deprecated, only use Eval
        void Sample(const shared_ptr<StageOCPControlSampler> &sampler, vector<double> &result);
        vector<double> Eval(const shared_ptr<StageOCPExprEvaluatorBase> &evaluator) const;
        void Eval(const shared_ptr<StageOCPExprEvaluatorBase> &evaluator, vector<double> &result) const;

    protected:
        vector<double> global_params;
        vector<double> stage_params;
        friend class StageOCPApplication;
    };

    class StageOCPApplication : public StageOCPApplicationAbstract
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

        private:
            const shared_ptr<OCPAdapter> adapter_;
        };

    public:
        shared_ptr<StageOCPExprEvaluatorFactory> GetExprEvaluator(const string &sampler_name);
        shared_ptr<StageOCPExprEvaluatorFactory> GetExprEvaluator(const shared_ptr<StageExpression> &expr)
        {
            return make_shared<StageOCPExprEvaluatorFactory>(expr, nu_, nx_, n_stage_params_, K_);
        }
        shared_ptr<AppParameterSetter> GetParameterSetter(const string &setter_name);
        void Build();
        int Optimize();
        const StageOCPSolution &LastStageOCPSolution();

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