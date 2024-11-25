#include "fatrop/ocp/OCPApplication.hpp"
#include "fatrop/solver/AlgBuilder.hpp"
#include "fatrop/ocp/OCPAdapter.hpp"
#include "fatrop/ocp/FatropOCP.hpp"
#include "fatrop/ocp/FatropOCPResto.hpp"
#include "fatrop/ocp/FatropOCPBuilder.hpp"
#include "fatrop/solver/FatropAlg.hpp"
#include "fatrop/ocp/OCPAbstract.hpp"
#include "fatrop/json/json.h"
#include "fatrop/auxiliary/Common.hpp"
#include "fatrop/solver/NLPL1.hpp"
using namespace fatrop;
using namespace std;
NLPApplication::NLPApplication() : fatropoptions_(make_shared<FatropOptions>()), dirty(true), options_registry(make_shared<FatropOptionsRegistry>(*fatropoptions_.get()))
{
    if (printer_ == nullptr)
    {
        printer_ = std::make_shared<FatropPrinter>();
    }
}

void NLPApplication::build(const shared_ptr<FatropNLP> &nlp, const shared_ptr<FatropNLP> &nlp_resto)
{
    // keep nlp around for getting nlpdims
    nlp_ = nlp;
    // check if prebuilt option "inequality_handling" is in prebuilt options
    if (fatropoptions_->inequality_handling.get() == "pd_ip")
    {
    }
    else if(fatropoptions_->inequality_handling.get() == "L1_pen")
    {
        nlp_ = std::make_shared<NLPL1>(nlp, fatropoptions_);
    }
    else
    {
        throw std::runtime_error("Unknown inequality handling method: " + fatropoptions_->inequality_handling.get());
    }
    std::shared_ptr<FatropData> fatropdata_resto;
    algbuilder = std::make_shared<AlgBuilder>();
    algbuilder -> set_printer(printer_);
    algbuilder -> build_fatrop_algorithm_objects(nlp_, nlp_resto, fatropoptions_, fatropdata_, fatropdata_resto, journaller_);
    printer_ -> print_level() = fatropoptions_->print_level.get();
    // fatropoptions_->register_option(IntegerOption::un_bounded("print_level", "prfatrop_fatrop_int level", &printer_->print_level(), 10));
    fatropalg_ = algbuilder -> build_algorithm();
    dirty = false;
}

fatrop_int NLPApplication::optimize() const
{
    assert(!dirty);
    // update options
    algbuilder -> update_options(*fatropoptions_);
    // solve optimization problem
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
    options_registry->set(option_name, value);
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
    shared_ptr<FatropNLP> nlp_resto(make_shared<FatropOCPResto>(FatropOCPBuilder(ocp_, fatropoptions_, printer_).build(ocp_), fatropoptions_));
    shared_ptr<FatropNLP> nlp(FatropOCPBuilder(ocp_, fatropoptions_, printer_).build(ocp_));
    NLPApplication::build(nlp, nlp_resto);
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