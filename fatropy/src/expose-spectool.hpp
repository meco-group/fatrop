#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "fatrop/ocp/StageOCPApplication.hpp"
#include "numpy.hpp"
#include "swig-bind.hpp"
#include "fatrop/spectool/spectool.hpp"

#include "fatrop/solver/AlgBuilder.hpp"
#include "fatrop/ocp/OCPAdapter.hpp"
#include "fatrop/ocp/FatropOCP.hpp"
#include "fatrop/ocp/FatropOCPBuilder.hpp"
#include "fatrop/ocp/StageOCP.hpp"
#include "fatrop/solver/FatropAlg.hpp"
#include "fatrop/ocp/OCPAbstract.hpp"
#include "fatrop/json/json.h"
#include "fatrop/auxiliary/Common.hpp"

namespace py = pybind11;
namespace fatropy
{
    struct casadi_Function_swig_wrap_type
    {
        typedef casadi::Function type;
        static constexpr char py_name[] = "Function";
        static constexpr char module[] = "casadi";
    };
    struct casadi_MX_swig_wrap_type
    {
        typedef casadi::MX type;
        static constexpr char py_name[] = "MX";
        static constexpr char module[] = "casadi";
    };
}

template <>
struct py::detail::type_caster<casadi::Function> : public py::detail::swig_type_caster<fatropy::casadi_Function_swig_wrap_type>
{
};
template <>
struct py::detail::type_caster<casadi::MX> : public py::detail::swig_type_caster<fatropy::casadi_MX_swig_wrap_type>
{
};
// template <>
// struct py::detail::type_caster<casadi::Dict>
// {

// };

// void print_function_nametest(const casadi::Function& f)
// {
//     std::cout << "Function pointer: ";
//     std::cout << &f << std::endl;
//     std::cout << "Function name: ";
//     std::cout << f.name() << std::endl;
//     std::cout << "number of inputs = " << f.n_in() << std::endl;
// };
// casadi::Function test_change_name(casadi::Function &f)
// {
//     return casadi::Function(f);
// }
// PYBIND11_MAKE_OPAQUE(casadi::Dict);
namespace fatropy
{
    struct ExposeSpectool
    {
        static void expose(py::module &m)
        {
            py::class_<casadi::GenericType>(m, "csGenericType")
                .def(py::init<bool>())
                .def(py::init<casadi_int>())
                .def(py::init<int>())
                .def(py::init<double>())
                .def(py::init<const std::string &>())
                .def(py::init<const std::vector<bool> &>())
                .def(py::init<const std::vector<casadi_int> &>())
                .def(py::init<const std::vector<int> &>())
                .def(py::init<const std::vector<std::vector<casadi_int>> &>())
                .def(py::init<const std::vector<double> &>())
                .def(py::init<const std::vector<std::vector<double>> &>())
                .def(py::init<const std::vector<std::string> &>())
                .def(py::init<const casadi::Function &>())
                .def(py::init<const std::vector<casadi::Function> &>())
                .def(py::init<const casadi::Dict &>())
                .def(py::init<void *>())
                .def("bool", &casadi::GenericType::operator bool)
                .def("int", &casadi::GenericType::operator casadi_int)
                .def("double", &casadi::GenericType::operator double)
                ;
            py::implicitly_convertible<bool, casadi::GenericType>();
            py::implicitly_convertible<casadi_int, casadi::GenericType>();
            py::implicitly_convertible<int, casadi::GenericType>();
            py::implicitly_convertible<double, casadi::GenericType>();
            py::implicitly_convertible<const std::string &, casadi::GenericType>();
            py::implicitly_convertible<const std::vector<bool> &, casadi::GenericType>();
            py::implicitly_convertible<const std::vector<casadi_int> &, casadi::GenericType>();
            py::implicitly_convertible<const std::vector<int> &, casadi::GenericType>();
            py::implicitly_convertible<const std::vector<std::vector<casadi_int>> &, casadi::GenericType>();
            py::implicitly_convertible<const std::vector<double> &, casadi::GenericType>();
            py::implicitly_convertible<const std::vector<std::vector<double>> &, casadi::GenericType>();
            py::implicitly_convertible<const std::vector<std::string> &, casadi::GenericType>();
            py::implicitly_convertible<const casadi::Function &, casadi::GenericType>();
            py::implicitly_convertible<const std::vector<casadi::Function> &, casadi::GenericType>();
            py::implicitly_convertible<const casadi::Dict &, casadi::GenericType>();
            py::bind_map<casadi::Dict>(m, "csDict");
            // this can maybe be used for nesting of dicts
            //             // ...
            //         .def("__init__", [](Bar &self, py::tuple t) {
            //             if (py::len(t) != 2)
            //                 throw std::runtime_error("oh no");
            //             new(&self) Bar(t[0].cast<int>(), t[1].cast<int>());
            //         })
            // // ...
            //     py::implicitly_convertible<py::tuple, Bar>();
            // py::implicitly_convertible<casadi::Dict, casadi::GenericType>();
            // py::implicitly_convertible<py::float_, casadi::MX>();
            // py::implicitly_convertible<double, >();
            py::class_<fatrop::spectool::Jacobian>(m, "Jacobian").def(py::init<casadi::MX, casadi::MX>());
            py::class_<fatrop::spectool::Hessian>(m, "Hessian").def(py::init<casadi::MX, casadi::MX, casadi::MX, casadi::MX>());

            py::bind_map<fatrop::spectool::uo_map_mx_mx>(m, "uo_map_mx_mx");
            py::class_<fatrop::spectool::IntegratorRk4>(m, "IntegratorRk4").def(py::init<const std::vector<std::pair<casadi::MX, casadi::MX>> &, const casadi::MX &>())
            .def("__call__",py::overload_cast<const casadi::MX&>(&fatrop::spectool::IntegratorRk4::operator()))
            .def("__call__",py::overload_cast<const std::vector<casadi::MX>&>(&fatrop::spectool::IntegratorRk4::operator()));
            py::enum_<fatrop::spectool::at>(m, "at")
                .value("t0", fatrop::spectool::at::t0)
                .value("mid", fatrop::spectool::at::mid)
                .value("tf", fatrop::spectool::at::tf)
                .export_values();
            py::class_<fatrop::spectool::Stage>(m, "Stage")
            .def("set_next", py::overload_cast<const casadi::MX &, const casadi::MX &, const fatrop::spectool::Jacobian &, const fatrop::spectool::Hessian &>(&fatrop::spectool::Stage::set_next), py::arg("x_next"), py::arg("expr"), py::arg("jacobian") = fatrop::spectool::Jacobian{}, py::arg("hessian") = fatrop::spectool::Hessian{})
            .def("add_objective", &fatrop::spectool::Stage::add_objective<fatrop::spectool::at>).def("add_objective", &fatrop::spectool::Stage::add_objective<fatrop::spectool::at, fatrop::spectool::at>).def("add_objective", &fatrop::spectool::Stage::add_objective<fatrop::spectool::at, fatrop::spectool::at, fatrop::spectool::at>).def("subject_to", &fatrop::spectool::Stage::subject_to<fatrop::spectool::at>).def("subject_to", &fatrop::spectool::Stage::subject_to<fatrop::spectool::at, fatrop::spectool::at>).def("subject_to", &fatrop::spectool::Stage::subject_to<fatrop::spectool::at, fatrop::spectool::at, fatrop::spectool::at>)
            .def("at_t0", py::overload_cast<const casadi::MX &>(&fatrop::spectool::Stage::at_t0, py::const_))
            .def("at_tf", py::overload_cast<const casadi::MX &>(&fatrop::spectool::Stage::at_tf, py::const_))
            .def("at_t0", py::overload_cast<>(&fatrop::spectool::Stage::at_t0, py::const_))
            .def("at_tf", py::overload_cast<>(&fatrop::spectool::Stage::at_tf, py::const_))
            .def("at_mid", &fatrop::spectool::Stage::at_mid)
            .def("sample", &fatrop::spectool::Stage::sample);


            // def("add_objective", &fatrop::spectool::Stage::add_objective<fatrop::spectool::at>);

            py::class_<fatrop::spectool::uStage>(m, "uStage")
                .def(py::init<>())
                .def(py::init<const int>())
                .def("add_objective", &fatrop::spectool::uStage::add_objective)
                .def("subject_to", py::overload_cast<const casadi::MX &, const fatrop::spectool::Jacobian &, const fatrop::spectool::Hessian &>(&fatrop::spectool::uStage::subject_to), py::arg("constraint"), py::arg("jacobian") = fatrop::spectool::Jacobian{}, py::arg("hessian") = fatrop::spectool::Hessian{})
                .def("set_next", py::overload_cast<const casadi::MX &, const casadi::MX &, const fatrop::spectool::Jacobian &, const fatrop::spectool::Hessian &>(&fatrop::spectool::uStage::set_next), py::arg("x_next"), py::arg("expr"), py::arg("jacobian") = fatrop::spectool::Jacobian{}, py::arg("hessian") = fatrop::spectool::Hessian{})
                // .def("set_next", py::overload_cast<const fatrop::spectool::uo_map_mx<casadi::MX>&>(&fatrop::spectool::Stage::set_next))
                .def("at_t0", &fatrop::spectool::uStage::at_t0)
                .def("at_tf", &fatrop::spectool::uStage::at_tf)
                .def("K", &fatrop::spectool::uStage::K)
                .def("dynamics", &fatrop::spectool::uStage::dynamics)
                .def("eval_at_control", &fatrop::spectool::uStage::eval_at_control)
                .def("sample", &fatrop::spectool::uStage::sample)
                .def("clone", &fatrop::spectool::uStage::clone)
                .def("register_state", &fatrop::spectool::uStage::register_state)
                .def("register_control", &fatrop::spectool::uStage::register_control)
                .def("register_hybrid", &fatrop::spectool::uStage::register_hybrid)
                .def("register_control_parameter", &fatrop::spectool::uStage::register_control_parameter)
                .def("register_global_parameter", &fatrop::spectool::uStage::register_global_parameter);
            py::class_<fatrop::spectool::SolverInterface, std::shared_ptr<fatrop::spectool::SolverInterface>>(m, "SolverInterface")
                .def("stats", &fatrop::spectool::SolverInterface::stats);
            // py::class_<std::shared_ptr<fatrop::spectool::SolverInterface>>(m, "SolverInterface_ptr")
            //     .def("stats", [](std::shared_ptr<fatrop::spectool::SolverInterface>& self){return self.get()->stats();});

            py::class_<fatrop::spectool::Ocp>(m, "Ocp")
                .def(py::init<>())
                .def("state", &fatrop::spectool::Ocp::state, py::arg("m") = 1, py::arg("n") = 1)
                .def("control", &fatrop::spectool::Ocp::control, py::arg("m") = 1, py::arg("n") = 1)
                .def("hybrid", &fatrop::spectool::Ocp::hybrid, py::arg("m") = 1, py::arg("n") = 1)
                .def("parameter", &fatrop::spectool::Ocp::parameter, py::arg("m") = 1, py::arg("n") = 1, py::arg("grid") = "global")
                .def("sample", &fatrop::spectool::Ocp::sample)
                .def("new_stage", &fatrop::spectool::Ocp::new_stage, py::arg("K") = 1)
                .def("new_ustage", &fatrop::spectool::Ocp::new_ustage, py::arg("K") = 1)
                .def("add_ustage", &fatrop::spectool::Ocp::add_ustage)
                .def("to_function", &fatrop::spectool::Ocp::to_function, py::arg("name"), py::arg("in"), py::arg("out"))
                .def("at_t0", py::overload_cast<const casadi::MX &>(&fatrop::spectool::Ocp::at_t0, py::const_))
                .def("at_tf", py::overload_cast<const casadi::MX &>(&fatrop::spectool::Ocp::at_tf, py::const_))
                .def("at_t0", py::overload_cast<>(&fatrop::spectool::Ocp::at_t0, py::const_))
                .def("at_tf", py::overload_cast<>(&fatrop::spectool::Ocp::at_tf, py::const_))
                .def("set_initial", &fatrop::spectool::Ocp::set_initial)
                .def("all_variables", &fatrop::spectool::Ocp::all_variables)
                .def("eval_at_initial", &fatrop::spectool::Ocp::eval_at_initial)
                .def("solver", &fatrop::spectool::Ocp::solver, py::arg("solver name"), py::arg("opts_function") = casadi::Dict(), py::arg("opts_fatrop") = casadi::Dict())
                .def("subject_to", &fatrop::spectool::Ocp::subject_to)
                .def("get_solver", &fatrop::spectool::Ocp::get_solver);
            // .def("state", &fatrop::spectool::Ocp::state, py::arg("m") = 1, py::arg("n") = 1);
        }
    };
}