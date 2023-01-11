#include "OCPBuilder.hpp"
using namespace fatrop;

OCPSolutionSampler::OCPSolutionSampler(int nu, int nx, int no_stage_params, int K, const shared_ptr<StageEvaluator> &eval, const shared_ptr<FatropData> &fatropdata, const shared_ptr<BFOCPAdapter> &ocp) : nu(nu),
                                                                                                                                                                                                            nx(nx),
                                                                                                                                                                                                            no_stage_params(no_stage_params),
                                                                                                                                                                                                            K_(K),
                                                                                                                                                                                                            eval_(eval),
                                                                                                                                                                                                            fatropdata_(fatropdata),
                                                                                                                                                                                                            ocp_(ocp)
{
}
int OCPSolutionSampler::Sample(vector<double> &sample)
{
    double *sol_p = ((VEC *)fatropdata_->x_curr)->pa;
    double *res_p = sample.data();
    double *global_params_p = ocp_->GetGlobalParams();
    double *stage_params_p = ocp_->GetStageParams();
    int size = eval_->Size();
    for (int k = 0; k < K_ - 1; k++)
    {
        eval_->Eval(sol_p + k * (nu + nx), sol_p + k * (nu + nx) + nu, global_params_p, stage_params_p + k * no_stage_params, res_p + k * size);
    };
    eval_->Eval(sol_p + (K_ - 2) * (nu + nx), sol_p + (K_ - 1) * (nu + nx), global_params_p, stage_params_p + (K_ - 1) * no_stage_params, res_p + (K_ - 1) * size);
    return 0;
}
vector<double> OCPSolutionSampler::Sample()
{
    vector<double> res(Size());
    Sample(res);
    return res;
}

OCPBuilder::OCPBuilder(const string &functions, const string &json_spec_file) : functions(functions), json_spec_file(json_spec_file)
{
}
shared_ptr<FatropApplication> OCPBuilder::Build()
{
    handle = make_shared<DLHandler>(functions);
    std::ifstream t(json_spec_file);
    std::stringstream buffer;
    buffer << t.rdbuf();
    json_spec = json::jobject::parse(buffer.str());
    K = json_spec["K"];
    nx = json_spec["nx"];
    nu = json_spec["nu"];
    const int ngI = json_spec["ngI"];
    const int ng = json_spec["ng"];
    const int ngF = json_spec["ngF"];
    const int ng_ineqI = json_spec["ng_ineqI"];
    const int ng_ineq = json_spec["ng_ineq"];
    const int ng_ineqF = json_spec["ng_ineqF"];
    no_stage_params = json_spec["n_stage_params"];
    no_global_params = json_spec["n_global_params"];
    EvalCasGen BAbtf(handle, "BAbt");
    EvalCasGen bkf(handle, "bk");
    EvalCasGen RSQrqtIf = GN ? EvalCasGen(handle, "RSQrqtIGN") : EvalCasGen(handle, "RSQrqtI");
    EvalCasGen rqIf(handle, "rqI");
    EvalCasGen RSQrqtf = GN ? EvalCasGen(handle, "RSQrqtGN") : EvalCasGen(handle, "RSQrqt");
    EvalCasGen rqf(handle, "rqk");
    EvalCasGen RSQrqtFf = GN ? EvalCasGen(handle, "RSQrqtFGN") : EvalCasGen(handle, "RSQrqtF");
    EvalCasGen rqFf(handle, "rqF");
    EvalCasGen GgtIf(handle, "GgtI");
    EvalCasGen gIf(handle, "gI");
    EvalCasGen Ggtf(handle, "Ggt");
    EvalCasGen gf(handle, "g");
    EvalCasGen GgtFf(handle, "GgtF");
    EvalCasGen gFf(handle, "gF");
    EvalCasGen LIf(handle, "LI");
    EvalCasGen Lkf(handle, "Lk");
    EvalCasGen LFf(handle, "LF");
    EvalCasGen GgineqItf(handle, "GgineqIt");
    EvalCasGen gineqIf(handle, "gineqI");
    EvalCasGen Ggineqtf(handle, "Ggineqt");
    EvalCasGen gineqf(handle, "gineq");
    EvalCasGen GgineqFtf(handle, "GgineqFt");
    EvalCasGen gineqFf(handle, "gineqF");
    ocptemplatebasic = make_shared<BFOCPBasic>(nu, nx, ngI, ng, ngF, ng_ineqI, ng_ineq, ng_ineqF, no_stage_params, no_global_params, K,
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
                                                                 Ggtf,
                                                                 gf,
                                                                 GgtFf,
                                                                 gFf,
                                                                 GgineqItf,
                                                                 gineqIf,
                                                                 Ggineqtf,
                                                                 gineqf,
                                                                 GgineqFtf,
                                                                 gineqFf,
                                                                 LIf,
                                                                 Lkf,
                                                                 LFf);
    ocptempladapteror = make_shared<BFOCPAdapter>(static_cast<shared_ptr<BFOCP>>(ocptemplatebasic));
    ocptempladapter = ocptempladapteror;
    ocptempladapter->SetParams(json_spec["stage_params"].get_number_array<double>("%lf"), json_spec["global_params"].get_number_array<double>("%lf"));
    shared_ptr<OCPLSRiccati> ocplsriccati1 = make_shared<OCPLSRiccati>(ocptempladapter->GetOCPDims());
    ocplsriccati = ocplsriccati1;
    params = make_shared<FatropParams>();
    ocpscaler = make_shared<OCPNoScaling>(params);
    shared_ptr<FatropOCP> fatropocp1 = make_shared<FatropOCP>(ocptempladapter, ocplsriccati, ocpscaler);
    fatropocp = fatropocp1;
    fatropdata = make_shared<FatropData>(fatropocp->GetNLPDims(), params);
    initial_u = json_spec["initial_u"].get_number_array<double>("%lf");
    initial_x = json_spec["initial_x"].get_number_array<double>("%lf");
    lowerI = json_spec["lowerI"].get_number_array<double>("%lf");
    upperI = json_spec["upperI"].get_number_array<double>("%lf");
    lower = json_spec["lower"].get_number_array<double>("%lf");
    upper = json_spec["upper"].get_number_array<double>("%lf");
    lowerF = json_spec["lowerF"].get_number_array<double>("%lf");
    upperF = json_spec["upperF"].get_number_array<double>("%lf");
    lower.insert(lower.begin(), lowerI.begin(), lowerI.end());
    upper.insert(upper.begin(), upperI.begin(), upperI.end());
    lower.insert(lower.end(), lowerF.begin(), lowerF.end());
    upper.insert(upper.end(), upperF.begin(), upperF.end());
    // std::vector<int> test = json_spec["states_offset"].as_object()["x1"].as_object().array(0).get_number_array<int>("%d");
    // auto test2 = json_spec["states_offset"][0];
    SetBounds();
    SetInitial();
    // vector<double> upper = vector<double>(lower.size(), INFINITY);
    filter = make_shared<Filter>(params->maxiter + 1);
    journaller = make_shared<Journaller>(params->maxiter + 1);
    linesearch = DDP ? make_shared<LineSearchDDP>(params, fatropocp, fatropdata, filter, journaller, ocplsriccati1, &(fatropocp1->ocpkktmemory_), ocptempladapter) : make_shared<BackTrackingLineSearch>(params, fatropocp, fatropdata, filter, journaller);
    fatropalg = make_shared<FatropAlg>(fatropocp, fatropdata, params, filter, linesearch, journaller);
    solver_built = true;
    vector<string> sampler_names = json_spec["samplers"];
    for (auto sampler_name : sampler_names)
    {
        sampler_map[sampler_name] = make_shared<OCPSolutionSampler>(GetSamplerCustom(sampler_name));
    }
    for (auto var_name: json_spec["states_offset"].as_object().keys())
    {
        sampler_map[string("state_") + var_name] = make_shared<OCPSolutionSampler>(GetSamplerState(var_name));
    }
    for (auto var_name: json_spec["controls_offset"].as_object().keys())
    {
        sampler_map[string("control_") + var_name] = make_shared<OCPSolutionSampler>(GetSamplerControl(var_name));
    }

    for (auto param_name: json_spec["global_params_offset"].as_object().keys())
    {
        parameter_setter_map[param_name] = make_shared<ParameterSetter>(GetParameterSetterGlobal(param_name));
    }
    for (auto param_name: json_spec["control_params_offset"].as_object().keys())
    {
        parameter_setter_map[param_name] = make_shared<ParameterSetter>(GetParameterSetterControl(param_name));
    }
    return fatropalg;
}
void OCPBuilder::SetBounds()
{
    // assert(solver_built);
    fatropdata->SetBounds(lower, upper);
}
void OCPBuilder::SetInitial()
{
    // assert(solver_built);
    // for(int k =0; k<K; k++)
    // {
    //     ocptemplatebasic->set_initial_uk(initial_u.data()+k*nu, k);
    //     ocptemplatebasic->set_initial_xk(initial_x.data()+k*nx, k);
    // }
    ocptempladapter->SetInitial(fatropdata, initial_u, initial_x);
}

int OCPBuilder::GetVariableMap(const string &variable_type, const string &variable_name, vector<int> &from, vector<int> &to)
{
    assert(solver_built);
    from = json_spec[variable_type].as_object()[variable_name].as_object().array(0).get_number_array<int>("%d");
    to = json_spec[variable_type].as_object()[variable_name].as_object().array(1).get_number_array<int>("%d");
    return 0;
}

int OCPBuilder::GetVariableMapState(const string &variable_name, vector<int> &from, vector<int> &to)
{
    assert(solver_built);
    return GetVariableMap("states_offset", variable_name, from, to);
}
int OCPBuilder::GetVariableMapControl(const string &variable_name, vector<int> &from, vector<int> &to)
{
    assert(solver_built);
    return GetVariableMap("controls_offset", variable_name, from, to);
}
int OCPBuilder::GetVariableMapControlParam(const string &variable_name, vector<int> &from, vector<int> &to)
{
    assert(solver_built);
    return GetVariableMap("control_params_offset", variable_name, from, to);
}
int OCPBuilder::GetVariableMapGlobalParam(const string &variable_name, vector<int> &from, vector<int> &to)
{
    assert(solver_built);
    return GetVariableMap("global_params_offset", variable_name, from, to);
}

OCPSolutionSampler OCPBuilder::GetSamplerState(const string &variable_name)
{
    assert(solver_built);
    vector<int> in;
    vector<int> out;
    GetVariableMapState(variable_name, in, out);
    return OCPSolutionSampler(nu, nx, no_stage_params, K, make_shared<IndexEvaluator>(false, in, out), fatropdata, ocptempladapteror);
}
OCPSolutionSampler OCPBuilder::GetSamplerControl(const string &variable_name)
{
    assert(solver_built);
    vector<int> in;
    vector<int> out;
    GetVariableMapControl(variable_name, in, out);
    return OCPSolutionSampler(nu, nx, no_stage_params, K, make_shared<IndexEvaluator>(true, in, out), fatropdata, ocptempladapteror);
}