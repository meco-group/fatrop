#ifndef OCPBUILDERALINCLUDED
#define OCPBUILDERALINCLUDED
#include "ocp/BFOCPBasic.hpp"
#include "ocp/BFOCPAL.hpp"
#include "ocp/BFOCPAdapter.hpp"
#include "ocp/BFOCPAdapterAL.hpp"
#include "ocp/OCPLSRiccati.hpp"
#include "ocp/OCPNoScaling.hpp"
#include "solver/FatropOptions.hpp"
#include "solver/Filter.hpp"
#include "ocp/FatropOCP.hpp"
#include "ocp/FatropOCPAL.hpp"
#include "solver/FatropAlg.hpp"
#include "alm_solver/FatropALMAlg.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include "json/json.h"
#include <sstream>
#include <templates/FatropNLPAL.hpp>
#include <templates/FatropApplication.hpp>
namespace fatrop
{
    class OCPBuilderAL
    {
    public:
        OCPBuilderAL(const string &functions, const string &json_spec_file);
        shared_ptr<FatropApplication> Build();
        const string functions;
        const string json_spec_file;
        bool GN = false;
        bool DDP = false;
        int K;
        shared_ptr<OCP> ocptempladapter;
        shared_ptr<BFOCPAdapterAL> ocptempladapterALor;
        shared_ptr<OCPAL> ocptempladapterAL;
        shared_ptr<OCPLinearSolver> ocplsriccati;
        shared_ptr<FatropOptions> params;
        shared_ptr<OCPScalingMethod> ocpscaler;
        shared_ptr<FatropNLPAL> fatropocpal;
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
        shared_ptr<FatropALMAlg> fatropalg;
        vector<double>& GlobalParams()
        {
            return ocptempladapterALor->globalparams;
        }
        private:
        void SetBounds();
        void SetInitial();
    };
}
#endif // OCPBUILDERINCLUDED