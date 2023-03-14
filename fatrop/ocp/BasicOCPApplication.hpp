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
        NLPApplication() : fatropparams_(make_shared<FatropParams>()), journaller_(make_shared<Journaller>(fatropparams_->maxiter + 1))
        {
        }

    protected:
        void Build(const shared_ptr<FatropNLP> &nlp)
        {
            // keep nlp around for getting nlpdims
            nlp_ = nlp;
            AlgBuilder algbuilder;
            algbuilder.BuildFatropAlgObjects(nlp, fatropparams_, fatropdata_, journaller_);
            fatropalg_ = algbuilder.BuildAlgorithm();
            dirty = false;
        }

    public:
        int Optimize()
        {
            assert(!dirty);
            int ret = fatropalg_->Optimize();
            return ret;
        }
        // TODO: make this protected and use last_solution instead and choose other name
        FatropVecBF &LastSolution()
        {
            assert(!dirty);
            return fatropdata_->x_curr;
        }
        FatropVecBF &InitialGuessPrimal()
        {
            assert(!dirty);
            return fatropdata_->x_initial;
        }
        FatropStats GetStats()
        {
            return fatropalg_->GetStats();
        }
        NLPDims GetNLPDims()
        {
            return nlp_->GetNLPDims();
        }

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
        OCPApplication(const shared_ptr<OCPAbstract> &ocp) : ocp_(ocp)
        {
        }

    protected:
        void Build()
        {
            // keep the adapter around for accessing the parameters for samplers and parameter setters
            adapter = make_shared<BFOCPAdapter>(ocp_);
            shared_ptr<FatropNLP> nlp(FatropOCPBuilder(ocp_, fatropparams_).Build(adapter));
            NLPApplication::Build(nlp);
            dirty = false;
        }

    public:
        vector<double> &GlobalParameters()
        {
            assert(!dirty);
            return adapter->GetGlobalParamsVec();
        }
        vector<double> &StageParameters()
        {
            assert(!dirty);
            return adapter->GetStageParamsVec();
        }
        void SetInitial(vector<double> &initial_u, vector<double> &initial_x)
        {
            assert(!dirty);
            adapter->SetInitial(fatropdata_, initial_u, initial_x);
        }
        OCPDims GetOCPDims()
        {
            return adapter->GetOCPDims();
        }

    private:
        const shared_ptr<OCPAbstract> ocp_;

    protected:
        shared_ptr<BFOCPAdapter> adapter;
    };

    class FatropSolution
    {
    public:
        void GetPrimalSolution(vector<double> &result)
        {
            result = sol_primal;
        }
        // defautl copy constructor
        FatropSolution(const FatropSolution &other) = default;

    protected:
        FatropSolution(){};
        void SetDims(const NLPDims &dims)
        {
            sol_primal.resize(dims.nvars);
        };
        void SetPrimalSolution(const FatropVecBF &sol)
        {
            sol.copyto(sol_primal);
        }

    protected:
        vector<double> sol_primal;
    };
    struct BasicOCPApplicationAbstract : public OCPApplication
    {
        BasicOCPApplicationAbstract(const shared_ptr<BasicOCP> &ocp) : OCPApplication(ocp)
        {
        }
    };

    struct BasicOCPSolution : public FatropSolution
    {
    public:
        BasicOCPSolution(const shared_ptr<BasicOCPApplicationAbstract> &app)
        {
            SetDims(app->GetOCPDims());
        }
        // void GetStates(vector<double> &result)
        // void GetControl(vector<double> &result)
    protected:
        BasicOCPSolution(){};
        int nx;
        int nu;
        int n_stage_params;
        int n_global_params;
        int K;
        void SetDims(const OCPDims &dims)
        {
            FatropSolution::SetDims(dims);
            nx = dims.nx.at(0);
            nu = dims.nu.at(0);
            n_stage_params = dims.n_stage_params.at(0);
            n_global_params = dims.n_global_params;
            K = dims.K;
            global_params.resize(n_global_params);
            stage_params.resize(n_stage_params);
        }
        void Set(const FatropVecBF &sol, const vector<double> &global_params, const vector<double> &stage_params)
        {
            FatropSolution::SetPrimalSolution(sol);
            this->global_params = global_params;
            this->stage_params = stage_params;
        }

    public:
        // todo make this deprecated, only use Eval
        void Sample(const shared_ptr<OCPControlSampler> &sampler, vector<double> &result)
        {
            sampler->Evaluate(sol_primal, global_params, stage_params, result);
        }
        vector<double> Eval(const shared_ptr<BasicOCPEvaluatorBase> &evaluator) const
        {
            return evaluator->Evaluate(sol_primal, global_params, stage_params);
        }
        void Eval(const shared_ptr<BasicOCPEvaluatorBase> &evaluator, vector<double> &result) const
        {
            evaluator->Evaluate(sol_primal, global_params, stage_params, result);
        }

    protected:
        vector<double> global_params;
        vector<double> stage_params;
        friend class BasicOCPApplication;
    };

    class BasicOCPApplication : public BasicOCPApplicationAbstract
    {
    public:
        BasicOCPApplication(const shared_ptr<BasicOCP> &ocp) : BasicOCPApplicationAbstract(ocp), nx_(ocp->nx_), nu_(ocp->nu_), K_(ocp->K_){};

    public:
        class AppParameterSetter : public ParameterSetter
        {
        public:
            AppParameterSetter(const shared_ptr<BFOCPAdapter> &adapter, const shared_ptr<ParameterSetter> &ps) : ParameterSetter(*ps), adapter_(adapter){};

        public:
            void SetValue(const double value[])
            {
                ParameterSetter::SetValue(adapter_->GetGlobalParamsVec(), adapter_->GetStageParamsVec(), value);
            };

        private:
            const shared_ptr<BFOCPAdapter> adapter_;
        };

    public:
        shared_ptr<BasicOCPEvaluatorFactory> GetEvaluator(const string &sampler_name)
        {
            return eval_factories[sampler_name];
        }
        shared_ptr<AppParameterSetter> GetParameterSetter(const string &setter_name)
        {
            return make_shared<AppParameterSetter>(adapter, param_setters[setter_name]);
        }
        void Build()
        {
            OCPApplication::Build();
            // allocate the last solution
            last_solution.SetDims(GetOCPDims());
            dirty = false;
        }
        int Optimize()
        {
            int ret = NLPApplication::Optimize();
            if (ret == 0)
            {
                last_solution.SetPrimalSolution(LastSolution());
            }
            return ret;
        }
        const BasicOCPSolution &LastBasicOCPSolution()
        {
            return last_solution;
        }

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
    { // auto eval = make_shared<EvalCasGen>(handle, "sampler_" + sampler_name);
      // result->samplers.insert(make_pair(sampler_name, make_shared<OCPSolutionSampler>(eval)));
    public:
        static shared_ptr<BasicOCPApplication> FromRockitInterface(const string &functions, const string &json_spec_file)
        {
            shared_ptr<DLHandler> handle = make_shared<DLHandler>(functions);
            std::ifstream t(json_spec_file);
            std::stringstream buffer;
            buffer << t.rdbuf();
            json::jobject json_spec = json::jobject::parse(buffer.str());
            auto ocptemplatebasic = BasicOCPBuilder::FromRockitInterface(handle, json_spec);
            // instantiate the BasicOCPApplication
            auto result = make_shared<BasicOCPApplication>(ocptemplatebasic);
            // add all samplers
            vector<string> sampler_names = json_spec["samplers"];
            const int nu = ocptemplatebasic->nu_;
            const int nx = ocptemplatebasic->nx_;
            const int no_stage_params = ocptemplatebasic->n_stage_params_;
            const int K = ocptemplatebasic->K_;
            for (auto sampler_name : sampler_names)
            {
                auto eval = make_shared<EvalCasGen>(handle, "sampler_" + sampler_name);
                result->eval_factories[sampler_name] = make_shared<BasicOCPEvaluatorFactory>(make_shared<EvalBaseSE>(eval), nu, nx, no_stage_params, K);
            }
            // add state samplers
            json::jobject states_offset = json_spec["states_offset"];
            vector<string> state_names = states_offset.keys();
            for (auto state_name : state_names)
            {
                vector<int> in = states_offset[state_name].as_object().array(0).get_number_array<int>("%d");
                vector<int> out = states_offset[state_name].as_object().array(1).get_number_array<int>("%d");
                result->eval_factories[string("state_") + state_name] = make_shared<BasicOCPEvaluatorFactory>(make_shared<IndexEvaluator>(false, in, out), nu, nx, no_stage_params, K);
            }
            // add control samplers
            json::jobject controls_offset = json_spec["controls_offset"];
            vector<string> control_names = controls_offset.keys();
            for (auto control_name : control_names)
            {
                vector<int> in = controls_offset[control_name].as_object().array(0).get_number_array<int>("%d");
                vector<int> out = controls_offset[control_name].as_object().array(1).get_number_array<int>("%d");
                result->eval_factories[string("control_") + control_name] = make_shared<BasicOCPEvaluatorFactory>(make_shared<IndexEvaluator>(true, in, out), nu, nx, no_stage_params, K);
            }
            // add all parameter setters
            json::jobject control_params_offset = json_spec["control_params_offset"];
            vector<string> control_params_names = control_params_offset.keys();
            for (auto control_params_name : control_params_names)
            {
                vector<int> in = control_params_offset[control_params_name].as_object().array(0).get_number_array<int>("%d");
                vector<int> out = control_params_offset[control_params_name].as_object().array(1).get_number_array<int>("%d");
                result->param_setters[control_params_name] = make_shared<ParameterSetter>(in, out, no_stage_params, in.size(), K, false);
            }
            json::jobject global_params_offset = json_spec["global_params_offset"];
            vector<string> global_params_names = global_params_offset.keys();
            for (auto global_params_name : global_params_names)
            {
                vector<int> in = global_params_offset[global_params_name].as_object().array(0).get_number_array<int>("%d");
                vector<int> out = global_params_offset[global_params_name].as_object().array(1).get_number_array<int>("%d");
                result->param_setters[global_params_name] = make_shared<ParameterSetter>(in, out, no_stage_params, in.size(), K, true);
            }
            result->Build();
            return result;
        }
    };
}
#endif // BASICOCPAPPLICATIONINCLUDED