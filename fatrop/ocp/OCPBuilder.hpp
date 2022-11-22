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
namespace fatrop
{
    class OCPSolutionSampler
    {
    public:
        OCPSolutionSampler(int nu, int nx, int K, int control, const vector<int> &offsets_in, const vector<int> &offsets_out, const shared_ptr<FatropData> &fatropdata);
        int Sample(vector<double> &sample);
    private:
        const int nu;
        const int nx;
        const int K;
        bool control;
        const vector<int> offsets_in;
        const vector<int> offsets_out;
        const int no_var;
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

    public:
        void SetBounds();
        void SetInitial();
        int GetVariableMap(const string &variable_type, const string &variable_name, vector<int> &from, vector<int> &to);
        int GetVariableMapStates(const string &variable_name, vector<int> &from, vector<int> &to);
        int GetVariableMapControls(const string &variable_name, vector<int> &from, vector<int> &to);
        int GetVariableMapControlParams(const string &variable_name, vector<int> &from, vector<int> &to);
        int GetVariableMapGlobalParams(const string &variable_name, vector<int> &from, vector<int> &to);
        OCPSolutionSampler GetSamplerState(const string &variable_name);
        OCPSolutionSampler GetSamplerControls(const string &variable_name);
    };
}
#endif // OCPBUILDERINCLUDED