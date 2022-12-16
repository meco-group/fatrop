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
namespace fatrop
{
    class StageEvaluator
    {
    public:
        virtual void Eval(double *u, double *x, double *global_params, double *stage_params, double *res) = 0;
        virtual int Size() = 0;
    };
    class IndexEvaluator : public StageEvaluator
    {
    public:
        IndexEvaluator(const bool control, const vector<int> offsets_in, const vector<int> offsets_out) : _no_var(offsets_in.size()),
                                                                                                          _offsets_in(offsets_in),
                                                                                                          _offsets_out(offsets_out),
                                                                                                          _control(control)
        {
        }
        void Eval(double *u, double *x, double *global_params, double *stage_params, double *res) override
        {
            if (_control)
            {
                for (int i = 0; i < _no_var; i++)
                {
                    res[_offsets_in.at(i)] = u[_offsets_out.at(i)];
                }
            }
            else
            {
                for (int i = 0; i < _no_var; i++)
                {
                    res[_offsets_in.at(i)] = x[_offsets_out.at(i)];
                }
            }
        };
        int Size() override
        {
            return _no_var;
        }

    private:
        const int _no_var;
        const vector<int> _offsets_in;
        const vector<int> _offsets_out;
        const bool _control;
    };
    class EvalBaseSE : public StageEvaluator
    {
    public:
        EvalBaseSE(const shared_ptr<EvalBase> &evalbase) : evalbase_(evalbase), size_(evalbase->out_m * evalbase->out_n)
        {
        }
        void Eval(double *u, double *x, double *global_params, double *stage_params, double *res) override
        {
            const double *arg[] = {u, x, stage_params, global_params};
            evalbase_->eval_array(arg, res);
        }
        int Size()
        {
            return size_;
        }

    private:
        shared_ptr<EvalBase> evalbase_;
        const int size_;
    };

    class OCPSolutionSampler
    {
    public:
        OCPSolutionSampler(int nu, int nx, int no_stage_params, int K, const shared_ptr<StageEvaluator> &eval, const shared_ptr<FatropData> &fatropdata);
        int Sample(vector<double> &sample);
        int Size()
        {
            return K*eval_->Size();
        }

    private:
        const int nu;
        const int nx;
        const int no_stage_params;
        const int K;
        shared_ptr<StageEvaluator> eval_;
        shared_ptr<FatropData> fatropdata_;
    };
    class OCPBuilder
    {
    public:
        OCPBuilder(const string &functions, const string &json_spec_file);
        shared_ptr<FatropApplication> Build();
        int K;
        int nu;
        int nx;
        const string functions;
        const string json_spec_file;
        bool solver_built = false;
        json::jobject json_spec;
        bool GN = false;
        bool DDP = false;
        shared_ptr<BFOCPAdapter> ocptempladapteror;
        shared_ptr<OCP> ocptempladapter;
        shared_ptr<OCPAL> ocptempladapterAL;
        shared_ptr<OCPLinearSolver> ocplsriccati;
        shared_ptr<FatropParams> params;
        shared_ptr<OCPScalingMethod> ocpscaler;
        shared_ptr<FatropNLP> fatropocp;
        shared_ptr<FatropData> fatropdata;
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

    public:
        void SetBounds();
        void SetInitial();
        int GetVariableMap(const string &variable_type, const string &variable_name, vector<int> &from, vector<int> &to);
        int GetVariableMapState(const string &variable_name, vector<int> &from, vector<int> &to);
        int GetVariableMapControl(const string &variable_name, vector<int> &from, vector<int> &to);
        int GetVariableMapControlParam(const string &variable_name, vector<int> &from, vector<int> &to);
        int GetVariableMapGlobalParam(const string &variable_name, vector<int> &from, vector<int> &to);
        OCPSolutionSampler GetSamplerState(const string &variable_name);
        OCPSolutionSampler GetSamplerControl(const string &variable_name);
        OCPSolutionSampler GetSamplerCustom(const string &sampler_name)
        {
            auto eval = make_shared<EvalCasGen>(handle, "sampler_" + sampler_name);
            return OCPSolutionSampler(nu, nx, 0, K, make_shared<EvalBaseSE>(eval), fatropdata);
        };
    };
}
#endif // OCPBUILDERINCLUDED