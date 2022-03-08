#ifndef OCPBUILDERINCLUDED
#define OCPBUILDERINCLUDED
#include "ocp/BFOCPBasic.hpp"
#include "ocp/BFOCPAdapter.hpp"
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
namespace fatrop
{
    class OCPBuilder
    {
    public:
        OCPBuilder(const string &functions, const string &json_spec_file)
        {
            std::ifstream t(json_spec_file);
            std::stringstream buffer;
            buffer << t.rdbuf();
            json::jobject json_spec = json::jobject::parse(buffer.str());
            K = json_spec["K"];
            const int nx = json_spec["nx"];
            const int nu = json_spec["nu"];
            const int ngI = json_spec["ngI"];
            const int ngF = json_spec["ngF"];
            const int ng_ineq = json_spec["ng_ineq"];
            const int ng_ineqF = json_spec["ng_ineqF"];
            const int n_stage_params = json_spec["n_stage_params"];
            const int n_global_params = json_spec["n_global_params"];
            RefCountPtr<DLHandler> handle = new DLHandler(functions);
            EvalCasGen BAbtf(handle, "BAbt");
            EvalCasGen bkf(handle, "bk");
            EvalCasGen RSQrqtIf(handle, "RSQrqtI");
            EvalCasGen rqIf(handle, "rqI");
            EvalCasGen RSQrqtf(handle, "RSQrqt");
            EvalCasGen rqf(handle, "rqk");
            EvalCasGen RSQrqtFf(handle, "RSQrqtF");
            EvalCasGen rqFf(handle, "rqF");
            EvalCasGen GgtIf(handle, "GgtI");
            EvalCasGen gIf(handle, "gI");
            EvalCasGen GgtFf(handle, "GgtF");
            EvalCasGen gFf(handle, "gF");
            EvalCasGen Lkf(handle, "Lk");
            EvalCasGen LFf(handle, "LF");
            EvalCasGen Ggineqtf(handle, "Ggineqt");
            EvalCasGen gineqf(handle, "gineq");
            EvalCasGen GgineqFtf(handle, "GgineqFt");
            EvalCasGen gineqFf(handle, "gineqF");
            RefCountPtr<BFOCP> ocptemplatebasic = new BFOCPBasic(nu, nx, ngI, ngF, ng_ineq, ng_ineqF, n_stage_params, n_global_params, K,
                                                                 BAbtf,
                                                                 bkf,
                                                                 RSQrqtIf,
                                                                 rqIf,
                                                                 RSQrqtf,
                                                                 rqf,
                                                                 RSQrqtFf,
                                                                 rqFf,
                                                                 GgtIf,
                                                                 gIf,
                                                                 GgtFf,
                                                                 gFf,
                                                                 Ggineqtf,
                                                                 gineqf,
                                                                 GgineqFtf,
                                                                 gineqFf,
                                                                 Lkf,
                                                                 LFf);
            ocptempladapter = new BFOCPAdapter(ocptemplatebasic);
            ocptempladapter->SetParams(json_spec["stage_params"].get_number_array<double>("%lf"), json_spec["global_params"].get_number_array<double>("%lf"));
            ocplsriccati = new OCPLSRiccati(ocptempladapter->GetOCPDims());
            params = new FatropParams();
            ocpscaler = new OCPNoScaling(params);
            fatropocp = new FatropOCP(ocptempladapter, ocplsriccati, ocpscaler);
            fatropdata = new FatropData(fatropocp->GetNLPDims(), params);
            initial_u = json_spec["initial_u"].get_number_array<double>("%lf");
            initial_x = json_spec["initial_x"].get_number_array<double>("%lf");
            lower = json_spec["lower"].get_number_array<double>("%lf");
            upper = json_spec["upper"].get_number_array<double>("%lf");
            lowerF = json_spec["lowerF"].get_number_array<double>("%lf");
            upperF = json_spec["upperF"].get_number_array<double>("%lf");
            lower.insert(lower.end(), lowerF.begin(), lowerF.end());
            upper.insert(upper.end(), upperF.begin(), upperF.end());
            SetBounds();
            SetInitial();
            // vector<double> upper = vector<double>(lower.size(), INFINITY);
            filter = new Filter(params->maxiter + 1);
            journaller = new Journaller(params->maxiter + 1);
            linesearch = new BackTrackingLineSearch(params, fatropocp, fatropdata, filter, journaller);
            fatropalg = new FatropAlg(fatropocp, fatropdata, params, filter, linesearch, journaller);
            // blasfeo_timer timer;
            // blasfeo_tic(&timer);
            // fatropalg->Optimize();
            // double el = blasfeo_toc(&timer);
            // cout << "el time " << el << endl;
        }
        void SetBounds()
        {
            ocptempladapter->SetInitial(K, fatropdata, initial_u, initial_x);
        }
        void SetInitial()
        {
            fatropdata->SetBounds(lower, upper);
        }
        int K;
        RefCountPtr<OCP> ocptempladapter;
        RefCountPtr<OCPLinearSolver> ocplsriccati;
        RefCountPtr<FatropParams> params;
        RefCountPtr<OCPScalingMethod> ocpscaler;
        RefCountPtr<FatropNLP> fatropocp;
        RefCountPtr<FatropData> fatropdata;
        vector<double> initial_u;
        vector<double> initial_x;
        vector<double> lower;
        vector<double> upper;
        vector<double> lowerF;
        vector<double> upperF;
        RefCountPtr<Filter> filter;
        RefCountPtr<Journaller> journaller;
        RefCountPtr<LineSearch> linesearch;
        RefCountPtr<FatropAlg> fatropalg;
    };
}
#endif // OCPBUILDERINCLUDED