#ifndef OCPBUILDERINCLUDED
#define OCPBUILDERINCLUDED
#include "ocp/BFOCPBasic.hpp"
#include "ocp/BFOCPAL.hpp"
#include "ocp/BFOCPAdapter.hpp"
#include "ocp/BFOCPAdapterAL.hpp"
#include "ocp/OCPLSRiccati.hpp"
#include "ocp/OCPNoScaling.hpp"
#include "solver/FatropParams.hpp"
#include "solver/Filter.hpp"
#include "ocp/FatropOCP.hpp"
#include "solver/FatropAlg.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include "json/json.h"
#include <sstream>
// #include <templates/FatropApplication.hpp>
#include <map>
#include "BasicOCPSamplers.hpp"
// #include <solver/AlgBuilder.hpp>
/*

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
THE OCPBUILDER CLASS IS DEPRECATED AS IT IS REPLACED BY THE BASICOCPAPPLICATION CLASS
It will be removed in the (near) future
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


*/
namespace fatrop
{
    class OCPSolutionSampler_old
    {
    public:
        OCPSolutionSampler_old(int nu, int nx, int no_stage_params, int K, const shared_ptr<StageEvaluator> &eval, const shared_ptr<FatropData> &fatropdata, const shared_ptr<BFOCPAdapter> &ocp);
        int Sample(vector<double> &sample);
        vector<double> Sample();
        int Size();
        int n_rows();
        int n_cols();
        int K();

    private:
        const int nu;
        const int nx;
        const int no_stage_params;
        const int K_;
        shared_ptr<StageEvaluator> eval_;
        shared_ptr<FatropData> fatropdata_;
        shared_ptr<BFOCPAdapter> ocp_;
    };

    class OCPBuilder
    {
    public:
        OCPBuilder(const string &functions, const string &json_spec_file);
        shared_ptr<FatropApplication> Build();
        int K;
        int nu;
        int nx;
        int no_global_params;
        int no_stage_params;
        const string functions;
        const string json_spec_file;
        bool solver_built = false;
        json::jobject json_spec;
        bool GN = false;
        bool DDP = false;
        shared_ptr<BasicOCP> ocptemplatebasic;
        shared_ptr<BFOCPAdapter> ocptempladapteror;
        shared_ptr<OCP> ocptempladapter;
        shared_ptr<OCPAL> ocptempladapterAL;
        shared_ptr<OCPLSRiccati> ocplsriccati1;
        shared_ptr<OCPLinearSolver> ocplsriccati;
        shared_ptr<FatropParams> params;
        shared_ptr<OCPScalingMethod> ocpscaler;
        shared_ptr<FatropNLP> fatropocp;
        shared_ptr<FatropData> fatropdata;
        shared_ptr<FatropOCP> fatropocp1;
        vector<double> initial_u;
        vector<double> initial_x;
        vector<double> lowerI;
        vector<double> upperI;
        vector<double> lower;
        vector<double> upper;
        vector<double> lowerF;
        vector<double> upperF;
        shared_ptr<Filter> filter;
        shared_ptr<Journaller> journaller;
        shared_ptr<LineSearch> linesearch;
        shared_ptr<FatropAlg> fatropalg;
        vector<double> &GlobalParams()
        {
            return ocptempladapteror->globalparams;
        }
        shared_ptr<DLHandler> handle;
        map<string, shared_ptr<OCPSolutionSampler_old>> sampler_map;
        map<string, shared_ptr<ParameterSetter>> parameter_setter_map;

    public:
        shared_ptr<OCPSolutionSampler_old> GetSampler(const string &sampler_name)
        {
            return sampler_map[sampler_name];
        }
        shared_ptr<ParameterSetter> GetParameterSetter(const string &parameter_setter_name)
        {
            return parameter_setter_map[parameter_setter_name];
        }

    private:
        void SetBounds();
        void SetInitial();
        int GetVariableMap(const string &variable_type, const string &variable_name, vector<int> &from, vector<int> &to);
        int GetVariableMapState(const string &variable_name, vector<int> &from, vector<int> &to);
        int GetVariableMapControl(const string &variable_name, vector<int> &from, vector<int> &to);
        int GetVariableMapControlParam(const string &variable_name, vector<int> &from, vector<int> &to);
        int GetVariableMapGlobalParam(const string &variable_name, vector<int> &from, vector<int> &to);
        ParameterSetter GetParameterSetterGlobal(const string &parameter_name)
        {
            vector<int> in;
            vector<int> out;
            GetVariableMapGlobalParam(parameter_name, in, out);
            return ParameterSetter(ocptempladapteror, in, out, no_stage_params, in.size(), K, true);
        }
        ParameterSetter GetParameterSetterControl(const string &parameter_name)
        {
            vector<int> in;
            vector<int> out;
            GetVariableMapControlParam(parameter_name, in, out);
            return ParameterSetter(ocptempladapteror, in, out, no_stage_params, in.size(), K, false);
        }
        OCPSolutionSampler_old GetSamplerState(const string &variable_name);
        OCPSolutionSampler_old GetSamplerControl(const string &variable_name);
        OCPSolutionSampler_old GetSamplerCustom(const string &sampler_name)
        {
            auto eval = make_shared<EvalCasGen>(handle, "sampler_" + sampler_name);
            return OCPSolutionSampler_old(nu, nx, no_stage_params, K, make_shared<EvalBaseSE>(eval), fatropdata, ocptempladapteror);
        };
    };
}
#endif // OCPBUILDERINCLUDED