#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
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
    struct casadi_dict_swig_wrap_type
    {
        typedef casadi::Dict type;
        static constexpr char py_name[] = "Dict";
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
template <>
struct py::detail::type_caster<casadi::Dict> : public py::detail::swig_type_caster<fatropy::casadi_dict_swig_wrap_type>
{
};

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
            py::class_<fatrop::spectrop::uo_map_mx<casadi::MX>>(m, "MXMXmap");
            py::class_<fatrop::spectrop::Stage>(m, "Stage")
            .def("add_objective", &fatrop::spectrop::Stage::add_objective)
            .def("subject_to", &fatrop::spectrop::Stage::subject_to)
            .def("set_next", py::overload_cast<const casadi::MX &, const casadi::MX& >(&fatrop::spectrop::Stage::set_next))
            // .def("set_next", py::overload_cast<const fatrop::spectrop::uo_map_mx<casadi::MX>&>(&fatrop::spectrop::Stage::set_next))
            .def("at_t0", &fatrop::spectrop::Stage::at_t0)
            .def("at_tf", &fatrop::spectrop::Stage::at_tf)
            .def("K", &fatrop::spectrop::Stage::K)
            .def("dynamics", &fatrop::spectrop::Stage::dynamics)
            .def("eval_at_control", &fatrop::spectrop::Stage::eval_at_control)
            .def("sample", &fatrop::spectrop::Stage::sample);

            py::class_<fatrop::spectrop::Ocp>(m, "Ocp").def(py::init<>())
            .def("state", &fatrop::spectrop::Ocp::state, py::arg("m") = 1, py::arg("n") = 1)
            .def("control", &fatrop::spectrop::Ocp::control, py::arg("m") = 1, py::arg("n") = 1)
            .def("automatic", &fatrop::spectrop::Ocp::automatic, py::arg("m") = 1, py::arg("n") = 1)
            .def("parameter", &fatrop::spectrop::Ocp::parameter, py::arg("m") = 1, py::arg("n") = 1, py::arg("grid") = "global")
            .def("sample", &fatrop::spectrop::Ocp::sample)
            .def("new_stage", &fatrop::spectrop::Ocp::new_stage, py::arg("K") = 1)
            .def("to_function", &fatrop::spectrop::Ocp::to_function, py::arg("from"), py::arg("to"), py::arg("opts") = casadi::Dict())
            .def("at_t0", &fatrop::spectrop::Ocp::at_t0)
            .def("at_tf", &fatrop::spectrop::Ocp::at_tf)
            .def("set_initial", &fatrop::spectrop::Ocp::set_initial);

            // .def("state", &fatrop::spectrop::Ocp::state, py::arg("m") = 1, py::arg("n") = 1);
        }
    };
}