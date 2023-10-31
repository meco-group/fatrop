#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "ocp/StageOCPApplication.hpp"
#include "numpy.hpp"
#include "swig-bind.hpp"
#include "spectrop.hpp"

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
namespace fatropy
{
    struct ExposeSpectrop
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
                .def(py::init<void *>());
            py::implicitly_convertible<bool, casadi::GenericType>();
            py::implicitly_convertible<casadi_int, casadi::GenericType>;
            py::implicitly_convertible<int, casadi::GenericType>;
            py::implicitly_convertible<double, casadi::GenericType>;
            py::implicitly_convertible<const std::string &, casadi::GenericType>;
            py::implicitly_convertible<const std::vector<bool> &, casadi::GenericType>;
            py::implicitly_convertible<const std::vector<casadi_int> &, casadi::GenericType>;
            py::implicitly_convertible<const std::vector<int> &, casadi::GenericType>;
            py::implicitly_convertible<const std::vector<std::vector<casadi_int>> &, casadi::GenericType>;
            py::implicitly_convertible<const std::vector<double> &, casadi::GenericType>;
            py::implicitly_convertible<const std::vector<std::vector<double>> &, casadi::GenericType>;
            py::implicitly_convertible<const std::vector<std::string> &, casadi::GenericType>;
            py::implicitly_convertible<const casadi::Function &, casadi::GenericType>;
            py::implicitly_convertible<const std::vector<casadi::Function> &, casadi::GenericType>;
            py::implicitly_convertible<const casadi::Dict &, casadi::GenericType>;
            // py::bind_map<casadi::Dict>(m, "cs::Dict");

            py::class_<fatrop::spectrop::uo_map_mx<casadi::MX>>(m, "MXMXmap");
            py::class_<fatrop::spectrop::uStage>(m, "Stage")
                .def("add_objective", &fatrop::spectrop::uStage::add_objective)
                .def("subject_to", &fatrop::spectrop::uStage::subject_to)
                .def("set_next", py::overload_cast<const casadi::MX &, const casadi::MX &>(&fatrop::spectrop::uStage::set_next))
                // .def("set_next", py::overload_cast<const fatrop::spectrop::uo_map_mx<casadi::MX>&>(&fatrop::spectrop::Stage::set_next))
                .def("at_t0", &fatrop::spectrop::uStage::at_t0)
                .def("at_tf", &fatrop::spectrop::uStage::at_tf)
                .def("K", &fatrop::spectrop::uStage::K)
                .def("dynamics", &fatrop::spectrop::uStage::dynamics)
                .def("eval_at_control", &fatrop::spectrop::uStage::eval_at_control)
                .def("sample", &fatrop::spectrop::uStage::sample);

            py::class_<fatrop::spectrop::Ocp>(m, "Ocp").def(py::init<>())
            .def("state", &fatrop::spectrop::Ocp::state, py::arg("m") = 1, py::arg("n") = 1)
            .def("control", &fatrop::spectrop::Ocp::control, py::arg("m") = 1, py::arg("n") = 1)
            .def("hybrid", &fatrop::spectrop::Ocp::hybrid, py::arg("m") = 1, py::arg("n") = 1)
            .def("parameter", &fatrop::spectrop::Ocp::parameter, py::arg("m") = 1, py::arg("n") = 1, py::arg("grid") = "global")
            .def("sample", &fatrop::spectrop::Ocp::sample)
            .def("new_stage", &fatrop::spectrop::Ocp::new_ustage, py::arg("K") = 1)
            .def("to_function", &fatrop::spectrop::Ocp::to_function, py::arg("in"), py::arg("out"), py::arg("opts") = casadi::Dict())
            .def("at_t0", &fatrop::spectrop::Ocp::at_t0)
            .def("at_tf", &fatrop::spectrop::Ocp::at_tf)
            .def("set_initial", &fatrop::spectrop::Ocp::set_initial);

            // .def("state", &fatrop::spectrop::Ocp::state, py::arg("m") = 1, py::arg("n") = 1);
        }
    };
}