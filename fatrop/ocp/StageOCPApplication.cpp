/*
 * Fatrop - A fast trajectory optimization solver
 *  Copyright (C) 2022 - 2024 Lander Vanroye, KU Leuven. All rights reserved.
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
#include "fatrop/ocp/StageOCPApplication.hpp"
#include "fatrop/solver/AlgBuilder.hpp"
#include "fatrop/ocp/OCPAdapter.hpp"
#include "fatrop/ocp/FatropOCP.hpp"
#include "fatrop/ocp/FatropOCPResto.hpp"
#include "fatrop/ocp/FatropOCPBuilder.hpp"
#include "fatrop/ocp/StageOCP.hpp"
#include "fatrop/solver/FatropAlg.hpp"
#include "fatrop/ocp/OCPAbstract.hpp"
#include "fatrop/json/json.h"
#include "fatrop/auxiliary/Common.hpp"
#include "fatrop/solver/NLPL1.hpp"
using namespace fatrop;
using namespace std;
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
StageOCPApplication::StageOCPApplication(const shared_ptr<StageOCP> &ocp) : OCPAbstractApplication(ocp), nx_(ocp->nx_), nu_(ocp->nu_), n_stage_params_(ocp->n_stage_params_), K_(ocp->K_){};
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
void StageOCPApplication::set_value(const std::string &setter_name, const std::vector<double> &value)
{
    get_parameter_setter(setter_name).set_value(value);
};
void StageOCPApplication::build()
{
    OCPAbstractApplication::build();
    // allocate the last solution
    last_solution_.set_dims(get_ocp_dims());
    dirty = false;
}
fatrop_int StageOCPApplication::optimize()
{
    fatrop_int ret = NLPApplication::optimize();
    last_solution_.set_parameters(global_parameters(), stage_parameters());
    last_solution_.set_solution(last_solution_primal(), last_solution_dual(), last_solution_zL(), last_solution_zU());
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