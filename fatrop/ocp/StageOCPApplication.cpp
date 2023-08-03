/*
 * Fatrop - A fast trajectory optimization solver
 * Copyright (C) 2022, 2023 Lander Vanroye, KU Leuven. All rights reserved.
 *
 * This file is part of Fatrop.
 *
 * Fatrop is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fatrop is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Fatrop.  If not, see <http://www.gnu.org/licenses/>. */
#include "ocp/StageOCPApplication.hpp"
using namespace fatrop;
using namespace std;
NLPApplication::NLPApplication() : fatropoptions_(make_shared<FatropOptions>()), dirty(true)
{
    if (printer_ == nullptr)
    {
        printer_ = std::make_shared<FatropPrinter>();
    }
}

void NLPApplication::build(const shared_ptr<FatropNLP> &nlp)
{
    // keep nlp around for getting nlpdims
    nlp_ = nlp;
    AlgBuilder algbuilder;
    algbuilder.set_printer(printer_);
    algbuilder.build_fatrop_algorithm_objects(nlp, fatropoptions_, fatropdata_, journaller_);
    fatropoptions_->register_option(IntegerOption::un_bounded("print_level", "prfatrop_fatrop_int level", &printer_->print_level(), 10));
    fatropalg_ = algbuilder.build_algorithm();
    dirty = false;
}

fatrop_int NLPApplication::optimize()
{
    assert(!dirty);
    fatrop_int ret = fatropalg_->optimize();
    return ret;
}
// TODO: make this protected and use last_solution instead and choose other name
const FatropVecBF &NLPApplication::last_solution_primal() const
{
    assert(!dirty);
    return fatropdata_->x_curr;
}
FatropVecBF &NLPApplication::initial_guess_primal() const
{
    assert(!dirty);
    return fatropdata_->x_initial;
}
FatropStats NLPApplication::get_stats() const
{
    return fatropalg_->get_stats();
}
NLPDims NLPApplication::get_nlp_dims()
{
    return nlp_->get_nlp_dims();
}
const FatropVecBF &NLPApplication::last_solution_dual() const
{
    return fatropdata_->lam_curr;
}
const FatropVecBF &NLPApplication::last_solution_zL() const
{
    return fatropdata_->zL_curr;
}
const FatropVecBF &NLPApplication::last_solution_zU() const
{
    return fatropdata_->zU_curr;
}
FatropVecBF &NLPApplication::initial_guess_dual() const
{
    return fatropdata_->lam_init;
}
FatropVecBF &NLPApplication::initial_guess_zL() const
{
    return fatropdata_->zL_init;
}
FatropVecBF &NLPApplication::initial_guess_zU() const
{
    return fatropdata_->zU_init;
}
template <typename T>
void NLPApplication::set_option(const string &option_name, T value)
{
    fatropoptions_->set(option_name, value);
}
template void NLPApplication::set_option<fatrop_int>(const string &, int);
template void NLPApplication::set_option<double>(const string &, double);
template void NLPApplication::set_option<bool>(const string &, bool);

void NLPApplication::set_initial(const FatropSolution &initial_guess)
{
    initial_guess_primal() = initial_guess.sol_primal_;
    initial_guess_dual() = initial_guess.sol_dual_;
    initial_guess_zL() = initial_guess.sol_zL_;
    initial_guess_zU() = initial_guess.sol_zU_;
}
const FatropOptions &NLPApplication::get_options() const
{
    return *fatropoptions_;
}

// TODO move this class to a separate file
OCPApplication::OCPApplication(const shared_ptr<OCPAbstract> &ocp) : ocp_(ocp)
{
}

void OCPApplication::build()
{
    // keep the adapter around for accessing the parameters for samplers and parameter setters
    adapter = make_shared<OCPAdapter>(ocp_, fatropoptions_);
    shared_ptr<FatropNLP> nlp(FatropOCPBuilder(ocp_, fatropoptions_, printer_).build(adapter));
    NLPApplication::build(nlp);
    dirty = false;
}

vector<double> &OCPApplication::global_parameters()
{
    assert(!dirty);
    return adapter->get_global_parameters_vec();
}
vector<double> &OCPApplication::stage_parameters()
{
    assert(!dirty);
    return adapter->get_stage_parameters_vec();
}
void OCPApplication::set_initial(vector<double> &initial_u, vector<double> &initial_x)
{
    assert(!dirty);
    adapter->set_initial_sol_guess(fatropdata_, initial_u, initial_x);
}
OCPDims OCPApplication::get_ocp_dims()
{
    return adapter->get_ocp_dims();
}

FatropSolution::FatropSolution(){};
void FatropSolution::set_dims(const NLPDims &dims)
{
    sol_primal_.resize(dims.nvars);
    sol_dual_.resize(dims.neqs);
    sol_zL_.resize(dims.nineqs);
    sol_zU_.resize(dims.nineqs);
};
void FatropSolution::set_solution(const FatropVecBF &sol_primal, const FatropVecBF &sol_dual, const FatropVecBF &sol_zL, const FatropVecBF &sol_zU)
{
    sol_primal.copyto(sol_primal_);
    sol_dual.copyto(sol_dual_);
    sol_zL.copyto(sol_zL_);
    sol_zU.copyto(sol_zU_);
};
void FatropSolution::set_primal_solution(const FatropVecBF &sol)
{
    sol.copyto(sol_primal_);
}
// StageOCPApplicationAbstract::StageOCPApplicationAbstract(const shared_ptr<StageOCP> &ocp) : OCPApplication(ocp)
// {
// }
StageOCPSolution::StageOCPSolution(const shared_ptr<OCP> &app)
{
    set_dims(app->get_ocp_dims());
}
StageOCPSolution::StageOCPSolution(){};
void StageOCPSolution::set_dims(const OCPDims &dims)
{
    FatropSolution::set_dims(dims);
    nx = dims.nx.get(0);
    nu = dims.nu.get(0);
    n_stage_params = dims.n_stage_params.get(0);
    n_global_params = dims.n_global_params;
    K = dims.K;
    global_params.resize(n_global_params);
    stage_params.resize(n_stage_params);
}
void StageOCPSolution::set_parameters(const vector<double> &global_params, const vector<double> &stage_params)
{
    this->global_params = global_params;
    this->stage_params = stage_params;
}

void StageOCPSolution::get_u(std::vector<double> &result) const
{
    result.resize(nu * (K - 1));
    for (fatrop_int k = 0; k < K - 1; k++)
    {
        fatrop_int offs = (nu + nx) * k;
        for (fatrop_int i = 0; i < nu; i++)
            result[nu * k + i] = sol_primal_[offs + i];
    }
}
void StageOCPSolution::get_x(std::vector<double> &result) const
{
    result.resize(nx * K);
    for (fatrop_int k = 0; k < K; k++)
    {
        fatrop_int offs = (k == K - 1) ? (nu + nx) * k : (nu + nx) * k + nu;
        for (fatrop_int i = 0; i < nx; i++)
            result[nx * k + i] = sol_primal_[offs + i];
    }
}
void StageOCPSolution::sample(const StageControlGridSampler &sampler, vector<double> &result) const
{
    sampler.evaluate(sol_primal_, global_params, stage_params, result);
}
vector<double> StageOCPSolution::evaluate(const StageExpressionEvaluatorBase &evaluator) const
{
    return evaluator.evaluate(sol_primal_, global_params, stage_params);
}
void StageOCPSolution::evaluate(const StageExpressionEvaluatorBase &evaluator, vector<double> &result) const
{
    evaluator.evaluate(sol_primal_, global_params, stage_params, result);
}
StageOCPApplication::StageOCPApplication(const shared_ptr<StageOCP> &ocp) : OCPApplication(ocp), nx_(ocp->nx_), nu_(ocp->nu_), n_stage_params_(ocp->n_stage_params_), K_(ocp->K_){};
void StageOCPApplication::set_initial_u(const std::vector<double> &initial_guess_u) const
{
    for (fatrop_int k = 0; k < K_ - 1; k++)
    {
        fatrop_int offs = (nu_ + nx_) * k;
        for (fatrop_int i = 0; i < nu_; i++)
            initial_guess_primal().at(offs + i) = initial_guess_u[nu_ * k + i];
    }
}
void StageOCPApplication::set_initial_x(const std::vector<double> &initial_guess_x) const
{
    for (fatrop_int k = 0; k < K_; k++)
    {
        fatrop_int offs = (k == K_ - 1) ? (nu_ + nx_) * k : (nu_ + nx_) * k + nu_;
        for (fatrop_int i = 0; i < nx_; i++)
            initial_guess_primal().at(offs + i) = initial_guess_x[nx_ * k + i];
    }
}

StageOCPApplication::AppParameterSetter::AppParameterSetter(const shared_ptr<OCPAdapter> &adapter, const shared_ptr<ParameterSetter> &ps) : ParameterSetter(*ps), adapter_(adapter){};
void StageOCPApplication::AppParameterSetter::set_value(const double value[])
{
    ParameterSetter::set_value(adapter_->get_global_parameters_vec(), adapter_->get_stage_parameters_vec(), value);
};

StageExpressionEvaluatorFactory StageOCPApplication::get_expression(const string &sampler_name)
{
    return get_expression_evaluator(stage_expressions[sampler_name]);
}
StageOCPApplication::AppParameterSetter StageOCPApplication::get_parameter_setter(const string &setter_name)
{
    return AppParameterSetter(adapter, param_setters[setter_name]);
}
void StageOCPApplication::build()
{
    OCPApplication::build();
    // allocate the last solution
    last_solution_.set_dims(get_ocp_dims());
    dirty = false;
}
fatrop_int StageOCPApplication::optimize()
{
    fatrop_int ret = NLPApplication::optimize();
    last_solution_.set_parameters(global_parameters(), stage_parameters());
    if (ret == 0)
    {
        last_solution_.set_solution(last_solution_primal(), last_solution_dual(), last_solution_zL(), last_solution_zU());
    }
    return ret;
}
const StageOCPSolution &StageOCPApplication::last_solution() const
{
    return last_solution_;
}

StageOCPApplication StageOCPApplicationFactory::from_rockit_interface(const string &functions, const string &json_spec_file)
{
    shared_ptr<DLHandler> handle = make_shared<DLHandler>(functions);
    std::ifstream t(json_spec_file);
    std::stringstream buffer;
    buffer << t.rdbuf();
    json::jobject json_spec = json::jobject::parse(buffer.str());
    auto stageocp = StageOCPBuilder::FromRockitInterface(handle, json_spec);
    // instantiate the BasicOCPApplication
    auto result = StageOCPApplication(stageocp);
    // add all samplers
    vector<string> sampler_names = json_spec["samplers"];
    // const fatrop_int nu = stageocp->nu_;
    // const fatrop_int nx = stageocp->nx_;
    const fatrop_int no_stage_params = stageocp->n_stage_params_;
    const fatrop_int K = stageocp->K_;
    for (auto sampler_name : sampler_names)
    {
        auto eval = make_shared<EvalCasGen>(handle, "sampler_" + sampler_name);
        result.stage_expressions[sampler_name] = make_shared<EvalBaseSE>(eval);
    }
    // add state samplers
    json::jobject states_offset = json_spec["states_offset"];
    vector<string> state_names = states_offset.keys();
    for (auto state_name : state_names)
    {
        vector<fatrop_int> in = states_offset[state_name].as_object().array(0).get_number_array<fatrop_int>("%d");
        vector<fatrop_int> out = states_offset[state_name].as_object().array(1).get_number_array<fatrop_int>("%d");
        result.stage_expressions[string("state_") + state_name] = make_shared<IndexEpression>(false, in, out);
    }
    // add control samplers
    json::jobject controls_offset = json_spec["controls_offset"];
    vector<string> control_names = controls_offset.keys();
    for (auto control_name : control_names)
    {
        vector<fatrop_int> in = controls_offset[control_name].as_object().array(0).get_number_array<fatrop_int>("%d");
        vector<fatrop_int> out = controls_offset[control_name].as_object().array(1).get_number_array<fatrop_int>("%d");
        result.stage_expressions[string("control_") + control_name] = make_shared<IndexEpression>(true, in, out);
    }
    // add all parameter setters
    json::jobject control_params_offset = json_spec["control_params_offset"];
    vector<string> control_params_names = control_params_offset.keys();
    for (auto control_params_name : control_params_names)
    {
        vector<fatrop_int> in = control_params_offset[control_params_name].as_object().array(0).get_number_array<fatrop_int>("%d");
        vector<fatrop_int> out = control_params_offset[control_params_name].as_object().array(1).get_number_array<fatrop_int>("%d");
        result.param_setters[control_params_name] = make_shared<ParameterSetter>(in, out, no_stage_params, in.size(), K, false);
    }
    json::jobject global_params_offset = json_spec["global_params_offset"];
    vector<string> global_params_names = global_params_offset.keys();
    for (auto global_params_name : global_params_names)
    {
        vector<fatrop_int> in = global_params_offset[global_params_name].as_object().array(0).get_number_array<fatrop_int>("%d");
        vector<fatrop_int> out = global_params_offset[global_params_name].as_object().array(1).get_number_array<fatrop_int>("%d");
        result.param_setters[global_params_name] = make_shared<ParameterSetter>(in, out, no_stage_params, in.size(), K, true);
    }
    result.build();
    return result;
}
void StageOCPApplication::AppParameterSetter::set_value(const initializer_list<double> il_)
{
    assert((int)il_.size() == _no_var);
    set_value(il_.begin());
}

const std::vector<std::string> StageOCPApplication::parameter_names() const
{
    std::vector<std::string> ret;
    for (auto &p : param_setters)
    {
        ret.push_back(p.first);
    }
    return ret;
}
const std::vector<std::string> StageOCPApplication::stage_expression_names() const
{
    std::vector<std::string> ret;
    for (auto &p : stage_expressions)
    {
        ret.push_back(p.first);
    }
    return ret;
}

StageExpressionEvaluatorFactory StageOCPApplication::get_expression_evaluator(const shared_ptr<StageExpression> &expr)
{
    return StageExpressionEvaluatorFactory(expr, nu_, nx_, n_stage_params_, K_);
}