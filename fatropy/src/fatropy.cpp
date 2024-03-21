#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "fatrop/ocp/StageOCPApplication.hpp"
#include "numpy.hpp"
#ifdef WITH_SPECTOOL
#include "swig-bind.hpp"
#include "specification.hpp"
#include "expose-spectool.hpp"
#endif

#include "fatrop/solver/AlgBuilder.hpp"
#include "fatrop/ocp/OCPAdapter.hpp"
#include "fatrop/ocp/FatropOCP.hpp"
#include "fatrop/ocp/FatropOCPBuilder.hpp"
#include "fatrop/ocp/StageOCP.hpp"
#include "fatrop/solver/FatropAlg.hpp"
#include "fatrop/ocp/OCPAbstract.hpp"
#include "fatrop/json/json.h"
#include "fatrop/auxiliary/Common.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
namespace fatropy
{
}

PYBIND11_MODULE(_fatropy, m)
{
    using namespace fatropy;
    // wrap print_function_name
    // py::class_<test_type>(m, "test_type").def(py::init<>());
    // m.def("print", &print);
    // m.def("print_function_name", [](const py::object &pyobj)
    //       {  print_function_nametest(*FromPySwig<casadi::Function>::convert(pyobj)); });
    // m.def("print_function_name", &print_function_nametest);
    // m.def("test_change_name", &test_change_name);
    py::class_<fatrop::FatropSolution>(m, "FatropSolution");
    py::class_<fatrop::OCPTimeStepSampler>(m, "OCPTimeStepSampler");
    py::class_<fatrop::StageControlGridSampler>(m, "StageControlGridSampler");
    py::class_<fatrop::FatropStats>(m, "FatropStats")
        .def_readonly("compute_sd_time", &fatrop::FatropStats::compute_sd_time)
        .def_readonly("duinf_time", &fatrop::FatropStats::duinf_time)
        .def_readonly("eval_hess_time", &fatrop::FatropStats::eval_hess_time)
        .def_readonly("eval_jac_time", &fatrop::FatropStats::eval_jac_time)
        .def_readonly("eval_cv_time", &fatrop::FatropStats::eval_cv_time)
        .def_readonly("eval_grad_time", &fatrop::FatropStats::eval_grad_time)
        .def_readonly("eval_obj_time", &fatrop::FatropStats::eval_obj_time)
        .def_readonly("initialization_time", &fatrop::FatropStats::initialization_time)
        .def_readonly("time_total", &fatrop::FatropStats::time_total)
        .def_readonly("eval_hess_count", &fatrop::FatropStats::eval_hess_count)
        .def_readonly("eval_jac_count", &fatrop::FatropStats::eval_jac_count)
        .def_readonly("eval_cv_count", &fatrop::FatropStats::eval_cv_count)
        .def_readonly("eval_grad_count", &fatrop::FatropStats::eval_grad_count)
        .def_readonly("eval_obj_count", &fatrop::FatropStats::eval_obj_count)
        .def_readonly("iterations_count", &fatrop::FatropStats::iterations_count);
    py::class_<fatrop::StageExpressionEvaluatorFactory>(m, "StageExpressionEvaluatorFactory")
        .def("at_t0", &fatrop::StageExpressionEvaluatorFactory::at_t0)
        .def("at_tf", &fatrop::StageExpressionEvaluatorFactory::at_tf)
        .def("at_control", &fatrop::StageExpressionEvaluatorFactory::at_control);
    py::class_<fatrop::StageOCPSolution, fatrop::FatropSolution>(m, "StageOCPSolution")
        .def(py::init<const std::shared_ptr<fatrop::OCP>>())
        .def_property_readonly("u", [](const fatrop::StageOCPSolution &sol)
                               {
            MatBind result;
            sol.get_u(result);
            return result.to_numpy(sol.get_nu(), sol.get_K()-1); })
        .def_property_readonly("x", [](const fatrop::StageOCPSolution &sol)
                               {
                MatBind result;
            sol.get_x(result);
            return result.to_numpy(sol.get_nx(), sol.get_K()); })
        .def("evaluate", [](const fatrop::StageOCPSolution &sol, const fatrop::OCPTimeStepSampler &t)
             {
                MatBind result;
                result.resize(t.size());
                sol.evaluate(t, result);
                return result.to_numpy(t.n_rows(), t.n_cols()); })
        .def("evaluate", [](const fatrop::StageOCPSolution &sol, const fatrop::StageControlGridSampler &t)
             {
                MatBind result;
                result.resize(t.size());
                sol.evaluate(t, result);
                if ((std::min)(t.n_rows(), t.n_cols()) == 1)
                    return result.to_numpy((std::max)(t.n_rows(), t.n_cols()), sol.get_K());
                else
                    return result.to_numpy(sol.get_K(), t.n_cols(), t.n_rows()); });

    py::class_<fatrop::StageOCPApplication>(m, "StageOCPApplication")
        .def(py::init<const std::shared_ptr<fatrop::StageOCP> &>())
        .def("optimize", &fatrop::StageOCPApplication::optimize)
        .def("set_option", &fatrop::StageOCPApplication::set_option<double>)
        .def("last_solution", &fatrop::StageOCPApplication::last_solution)
        .def("get_expression", &fatrop::StageOCPApplication::get_expression)
        // .def("set_initial", [](fatrop::NLPApplication& app,  const fatrop::FatropSolution &initial_guess){return app.set_initial(initial_guess);}, py::const_)
        .def("set_initial", py::overload_cast<const fatrop::FatropSolution &>(&fatrop::NLPApplication::set_initial, py::const_))
        .def("set_initial_x", [](fatrop::StageOCPApplication &app, const py::array_t<double> &x)
             { app.set_initial_x(MatBind(x)); })
        .def("set_initial_u", [](fatrop::StageOCPApplication &app, const py::array_t<double> &u)
             { app.set_initial_u(MatBind(u)); })
        .def("set_value",
             [](fatrop::StageOCPApplication &self, const std::string &name, const py::array_t<double> &value)
             {
                 self.set_value(name, MatBind(value));
             })
        .def("set_params",
             [](fatrop::StageOCPApplication &self, const py::array_t<double> &stage_params, const py::array_t<double> &global_params)
             {
                 self.set_params(MatBind(global_params), MatBind(stage_params));
             })
        .def("get_stats", &fatrop::StageOCPApplication::get_stats);

    py::class_<fatrop::StageOCPApplicationFactory>(m, "StageOCPApplicationFactory")
        .def(py::init<>())
        .def("from_rockit_interface", &fatrop::StageOCPApplicationFactory::from_rockit_interface);

#ifdef WITH_SPECTOOL
    auto module_spectool = m.def_submodule("spectool", "specification tool");
    ExposeSpectool::expose(module_spectool);
#endif

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
