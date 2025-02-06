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
#include "fatrop/ocp/StageOCPApplication.hpp"
#include "fatrop/solver/AlgBuilder.hpp"
#include "fatrop/ocp/OCPAdapter.hpp"
#include "fatrop/ocp/FatropOCP.hpp"
#include "fatrop/ocp/FatropOCPBuilder.hpp"
#include "fatrop/ocp/StageOCP.hpp"
#include "fatrop/solver/FatropAlg.hpp"
#include "fatrop/ocp/OCPAbstract.hpp"
#include "fatrop/auxiliary/Common.hpp"
#include "fatrop/solver/NLPL1.hpp"
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
    // check if prebuilt option "inequality_handling" is in prebuilt options
    if (fatropoptions_->prebuilt_string.find("inequality_handling") == fatropoptions_->prebuilt_string.end())
    {
        fatropoptions_ ->prebuilt_set<std::string>("inequality_handling", "pd_ip");
    }
    if(fatropoptions_->prebuilt_string["inequality_handling"] == "L1_pen")
    {
        nlp_ = std::make_shared<NLPL1>(nlp, fatropoptions_);
    }
    AlgBuilder algbuilder;
    algbuilder.set_printer(printer_);
    algbuilder.build_fatrop_algorithm_objects(nlp_, fatropoptions_, fatropdata_, journaller_);
    fatropoptions_->register_option(IntegerOption::un_bounded("print_level", "prfatrop_fatrop_int level", &printer_->print_level(), 10));
    fatropalg_ = algbuilder.build_algorithm();
    dirty = false;
}

fatrop_int NLPApplication::optimize() const
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
    // check if application is built
    if(dirty)
    fatropoptions_->prebuilt_set(option_name, value);
    else
    fatropoptions_->set(option_name, value);
}
template void NLPApplication::set_option<fatrop_int>(const string &, int);
template void NLPApplication::set_option<double>(const string &, double);
template void NLPApplication::set_option<bool>(const string &, bool);
template void NLPApplication::set_option<string>(const string &, string);

void NLPApplication::set_initial(const FatropSolution &initial_guess) const
{
    initial_guess_primal() = initial_guess.sol_primal_;
    initial_guess_dual() = initial_guess.sol_dual_;
    initial_guess_zL() = initial_guess.sol_zL_;
    initial_guess_zU() = initial_guess.sol_zU_;
}
void NLPApplication::set_initial(const std::vector<double> &initial_guess_primal_) const
{
    initial_guess_primal() = initial_guess_primal_;
}
const FatropOptions &NLPApplication::get_options() const
{
    return *fatropoptions_;
}

// TODO move this class to a separate file
OCPAbstractApplication::OCPAbstractApplication(const shared_ptr<OCPAbstract> &ocp) 
{
    adapter = make_shared<OCPAdapter>(ocp, fatropoptions_);
    ocp_ = adapter;
}

void OCPApplication::build()
{
    // keep the adapter around for accessing the parameters for samplers and parameter setters
    shared_ptr<FatropNLP> nlp(FatropOCPBuilder(ocp_, fatropoptions_, printer_).build(ocp_));
    NLPApplication::build(nlp);
    dirty = false;
}
void OCPAbstractApplication::set_params(const std::vector<double> &global_params, const std::vector<double> &stage_params)
{
    adapter->set_parameters(stage_params, global_params);
}

vector<double> &OCPAbstractApplication::global_parameters()
{
    assert(!dirty);
    return adapter->get_global_parameters_vec();
}
vector<double> &OCPAbstractApplication::stage_parameters()
{
    assert(!dirty);
    return adapter->get_stage_parameters_vec();
}


OCPApplication::OCPApplication(const std::shared_ptr<OCP> &ocp):ocp_(ocp)
{

}
OCPApplication::OCPApplication():ocp_(nullptr)
{

}
void OCPApplication::set_initial(vector<double> &initial_u, vector<double> &initial_x)
{
    assert(!dirty);
    ocp_->set_initial_sol_guess(fatropdata_, initial_u, initial_x);
}
OCPDims OCPApplication::get_ocp_dims()
{
    return ocp_->get_ocp_dims();
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
    sol_primal.block(0, sol_primal_.size()).copyto(sol_primal_);
    sol_dual.block(0, sol_dual_.size()).copyto(sol_dual_);
    sol_zL.block(0, sol_zL_.size()).copyto(sol_zL_);
    sol_zU.block(0,sol_zU_.size()).copyto(sol_zU_);
};
void FatropSolution::set_primal_solution(const FatropVecBF &sol)
{
    sol.block(0, sol_primal_.size()).copyto(sol_primal_);
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
