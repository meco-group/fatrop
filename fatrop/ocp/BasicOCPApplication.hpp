#ifndef BASICOCPAPPLICATIONINCLUDED
#define BASICOCPAPPLICATIONINCLUDED
#include "ocp/BFOCPBasic.hpp"
#include "solver/AlgBuilder.hpp"
#include "ocp/BFOCPAdapter.hpp"
#include "solver/FatropAlg.hpp"
#include "ocp/FatropOCP.hpp"
#include "ocp/OCPLSRiccati.hpp"
#include "ocp/OCPNoScaling.hpp"
#include "ocp/FatropOCPBuilder.hpp"
#include "ocp/OCPAbstract.hpp"
#include "ocp/BasicOCPSamplers.hpp"
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
        void Build(const shared_ptr<FatropNLP> &nlp)
        {
            AlgBuilder algbuilder;
            algbuilder.BuildFatropAlgObjects(nlp, fatropparams_, fatropdata_, journaller_);
            fatropalg_ = algbuilder.BuildAlgorithm();
            dirty = false;
        }
        int Optimize()
        {
            assert(!dirty);
            return fatropalg_->Optimize();
        }
        FatropVecBF &LastSolution()
        {
            assert(!dirty);
            return fatropdata_->x_curr;
        }
        bool dirty = true;
        const shared_ptr<OCPAbstract> ocp_;
        shared_ptr<FatropParams> fatropparams_;
        shared_ptr<Journaller> journaller_;
        shared_ptr<FatropData> fatropdata_;
        shared_ptr<FatropAlg> fatropalg_;
    };

    // TODO move this class to a separate file
    class OCPApplication : public NLPApplication
    {
    public:
        OCPApplication(const shared_ptr<OCPAbstract> &ocp) : ocp_(ocp)
        {
        }
        void Build()
        {
            // keep the adapter around for accessing the parameters for samplers and parameter setters
            adapter = make_shared<BFOCPAdapter>(ocp_);
            shared_ptr<FatropNLP> nlp(FatropOCPBuilder(ocp_, fatropparams_).Build(adapter));
            NLPApplication::Build(nlp);
            dirty = false;
        }
        int Optimize()
        {
            assert(!dirty);
            return NLPApplication::Optimize();
        }
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

    private:
        bool dirty = true;
        const shared_ptr<OCPAbstract> ocp_;
        shared_ptr<BFOCPAdapter> adapter;
    };

    class BasicOCPApplication : public OCPApplication
    {
    public:
        BasicOCPApplication(const shared_ptr<BasicOCP> &ocp) : OCPApplication(ocp){};
        shared_ptr<OCPSolutionSampler> GetSampler(const string &sampler_name)
        {
            return samplers[sampler_name];
        }
        map<string, shared_ptr<OCPSolutionSampler>> samplers;
    };

    class BasicOCPApplicationBuilder
    {                // auto eval = make_shared<EvalCasGen>(handle, "sampler_" + sampler_name);
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
            cout << "Found " << sampler_names.size() << " samplers" << endl;
            const int nu = ocptemplatebasic -> nu_;
            const int nx = ocptemplatebasic -> nx_;
            const int no_stage_params = ocptemplatebasic -> n_stage_params_;
            const int K = ocptemplatebasic->K_;
            for (auto sampler_name : sampler_names)
            {
                cout << sampler_name << endl;
                auto eval = make_shared<EvalCasGen>(handle, "sampler_" + sampler_name);
                result->samplers.insert(make_pair(sampler_name, make_shared<OCPSolutionSampler>(nu, nx, no_stage_params, K, make_shared<EvalBaseSE>(eval))));
            }
            // string sampler_name = "x";
            // auto eval = make_shared<EvalCasGen>(handle, "sampler_" + sampler_name);
            // result -> samplers.insert(make_pair(sampler_name, make_shared<OCPSolutionSampler>(eval)));
            // add all parameter setters
            return result;
        }
    };
}
#endif // BASICOCPAPPLICATIONINCLUDED