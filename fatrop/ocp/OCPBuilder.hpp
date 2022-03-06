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
            const int K = json_spec["K"];
            const int nx = json_spec["nx"];
            const int nu = json_spec["nu"];
            const int ngI = json_spec["ngI"];
            const int ngF = json_spec["ngF"];
            const int ng_ineq = json_spec["ng_ineq"];
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
            RefCountPtr<BFOCP> ocptemplatebasic = new BFOCPBasic(nu, nx, ngI, ngF, ng_ineq, n_stage_params, n_global_params, K,
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
                                                                 Lkf,
                                                                 LFf);
            RefCountPtr<OCP> ocptempladapter = new BFOCPAdapter(ocptemplatebasic);
            ocptempladapter->SetParams(json_spec["stage_params"].get_number_array<double>("%lf"), json_spec["global_params"].get_number_array<double>("%lf"));
            RefCountPtr<OCPLinearSolver> ocplsriccati = new OCPLSRiccati(ocptempladapter->GetOCPDims());
            RefCountPtr<FatropParams> params = new FatropParams();
            RefCountPtr<OCPScalingMethod> ocpscaler = new OCPNoScaling(params);
            RefCountPtr<FatropNLP> fatropocp = new FatropOCP(ocptempladapter, ocplsriccati, ocpscaler);
            RefCountPtr<FatropData> fatropdata = new FatropData(fatropocp->GetNLPDims(), params);
            vector<double> initial_u = json_spec["initial_u"].get_number_array<double>("%lf");
            vector<double> initial_x = json_spec["initial_x"].get_number_array<double>("%lf");
            vector<double> lower = json_spec["lower"].get_number_array<double>("%lf");
            // vector<double> upper = json_spec["upper"].get_number_array<double>("%lf");
            vector<double> upper = vector<double>(lower.size(), INFINITY);
            ocptempladapter->SetInitial(K, fatropdata, initial_u, initial_x);
            fatropdata ->SetBounds(lower, upper);
            RefCountPtr<Filter> filter(new Filter(params->maxiter + 1));
            RefCountPtr<Journaller> journaller(new Journaller(params->maxiter + 1));
            RefCountPtr<LineSearch> linesearch = new BackTrackingLineSearch(params, fatropocp, fatropdata, filter, journaller);
            RefCountPtr<FatropAlg> fatropalg = new FatropAlg(fatropocp, fatropdata, params, filter, linesearch, journaller);
            blasfeo_timer timer;
            blasfeo_tic(&timer);
            fatropalg->Optimize();
            double el = blasfeo_toc(&timer);
            cout << "el time " << el << endl;
        }
    };
}
#endif // OCPBUILDERINCLUDED