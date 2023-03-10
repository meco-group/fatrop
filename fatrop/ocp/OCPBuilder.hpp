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
#include <templates/FatropApplication.hpp>
#include <map>
// #include <solver/AlgBuilder.hpp>
namespace fatrop
{
    class StageEvaluator
    {
    public:
        virtual void Eval(double *u, double *x, double *global_params, double *stage_params, double *res) = 0;
        virtual int n_rows() = 0;
        virtual int n_cols() = 0;
        int Size();
    };
    class IndexEvaluator : public StageEvaluator
    {
    public:
        IndexEvaluator(const bool control, const vector<int> offsets_in, const vector<int> offsets_out);
        void Eval(double *u, double *x, double *global_params, double *stage_params, double *res)override;
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
        void Eval(double *u, double *x, double *global_params, double *stage_params, double *res) override;
        int n_rows() override;
        int n_cols() override;

    private:
        shared_ptr<EvalBase> evalbase_;
        const int n_rows_;
        const int n_cols_;
    };

    class OCPSolutionSampler
    {
    public:
        OCPSolutionSampler(int nu, int nx, int no_stage_params, int K, const shared_ptr<StageEvaluator> &eval, const shared_ptr<FatropData> &fatropdata, const shared_ptr<BFOCPAdapter> &ocp);
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

    class ParameterSetter
    {
    public:
        ParameterSetter(const shared_ptr<BFOCPAdapter> &ocp, const vector<int> &offsets_in, const vector<int> &offsets_out, const int no_stage_params, const int no_var, const int K, const bool global);
        void SetValue(const double value[]);
        void SetValue(const initializer_list<double> il_);

    private:
        shared_ptr<BFOCPAdapter> ocp_;
        const vector<int> _offsets_in;
        const vector<int> _offsets_out;
        const int no_stage_params;
        const int _no_var;
        const int K;
        const bool _global;
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
        map<string, shared_ptr<OCPSolutionSampler>> sampler_map;
        map<string, shared_ptr<ParameterSetter>> parameter_setter_map;

    public:
        shared_ptr<OCPSolutionSampler> GetSampler(const string &sampler_name)
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
        OCPSolutionSampler GetSamplerState(const string &variable_name);
        OCPSolutionSampler GetSamplerControl(const string &variable_name);
        OCPSolutionSampler GetSamplerCustom(const string &sampler_name)
        {
            auto eval = make_shared<EvalCasGen>(handle, "sampler_" + sampler_name);
            return OCPSolutionSampler(nu, nx, no_stage_params, K, make_shared<EvalBaseSE>(eval), fatropdata, ocptempladapteror);
        };
    };
}
#endif // OCPBUILDERINCLUDED