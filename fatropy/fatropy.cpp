#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "ocp/StageOCPApplication.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j)
{
    return i + j;
}

namespace py = pybind11;
namespace fatropy
{
    struct MatBind : public std::vector<double>
    {
        using std::vector<double>::vector;
        py::array_t<double> to_numpy()
        {
            auto result = py::array(py::buffer_info(
                nullptr,                              /* Pointer to data (nullptr -> ask NumPy to allocate!) */
                sizeof(double),                       /* Size of one item */
                py::format_descriptor<double>::value, /* Buffer format */
                1,                                    /* How many dimensions? */
                {this->size()},                       /* Number of elements for each dimension */
                {sizeof(double)}                      /* Strides for each dimension */
                ));
            std::copy(this->begin(), this->end(), reinterpret_cast<double *>(result.request().ptr));
            return result;
        }
        // 2d array to_numpy()
        py::array_t<double> to_numpy(const int m, const int n)
        {
            // result is in Fortran format
            auto result = py::array(py::buffer_info(
                nullptr,                              /* Pointer to data (nullptr -> ask NumPy to allocate!) */
                sizeof(double),                       /* Size of one item */
                py::format_descriptor<double>::value, /* Buffer format */
                2,                                    /* How many dimensions? */
                {m, n},                               /* Number of elements for each dimension */
                {sizeof(double), sizeof(double) * m}  /* Strides for each dimension */
                ));
            std::copy(this->begin(), this->end(), reinterpret_cast<double *>(result.request().ptr));
            return result;
        }
        py::array_t<double> to_numpy(const int l, const int m, const int n)
        {
            // l is the number of 2d matrices
            auto result = py::array(py::buffer_info(
                nullptr,                                                     /* Pointer to data (nullptr -> ask NumPy to allocate!) */
                sizeof(double),                                              /* Size of one item */
                py::format_descriptor<double>::value,                        /* Buffer format */
                3,                                                           /* How many dimensions? */
                {l, m, n},                                                   /* Number of elements for each dimension */
                {sizeof(double) * m * n, sizeof(double), sizeof(double) * m} /* Strides for each dimension */
                ));
            std::copy(this->begin(), this->end(), reinterpret_cast<double *>(result.request().ptr));
            return result;
        }
    };
}

PYBIND11_MODULE(fatropy, m)
{
    using namespace fatropy;
    py::class_<fatrop::FatropSolution>(m, "FatropSolution");
    py::class_<fatrop::OCPTimeStepSampler>(m, "OCPTimeStepSampler");
    py::class_<fatrop::StageControlGridSampler>(m, "StageControlGridSampler");
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
                if (std::min(t.n_rows(), t.n_cols()) == 1)
                    return result.to_numpy(std::max(t.n_rows(), t.n_cols()), sol.get_K());
                else
                    return result.to_numpy(sol.get_K(), t.n_cols(), t.n_rows()); });


    py::class_<fatrop::StageOCPApplication>(m, "StageOCPApplication")
        .def(py::init<const std::shared_ptr<fatrop::StageOCP> &>())
        .def("optimize", &fatrop::StageOCPApplication::optimize)
        .def("set_option", &fatrop::StageOCPApplication::set_option<double>)
        .def("last_stageocp_solution", &fatrop::StageOCPApplication::last_stageocp_solution)
        .def("get_expression", &fatrop::StageOCPApplication::get_expression)
        .def("set_initial", &fatrop::NLPApplication::set_initial)
        // .def("set_initial_x", [](fatrop::StageOCPApplication &app, const py::array_t<double> &x)
        //      {
        //         auto buf = x.request();
        //         if (buf.ndim != 2)
        //             throw std::runtime_error("Number of dimensions must be two");
        //         if (buf.shape[0] != app.get_nx())
        //             throw std::runtime_error("Number of rows must be equal to nx");
        //         if (buf.shape[1] != app.get_K())
        //             throw std::runtime_error("Number of columns must be equal to K");
        //         app.set_initial_x(reinterpret_cast<double *>(buf.ptr));
        //      });
        ;
        
    py::class_<fatrop::StageOCPApplicationFactory>(m, "StageOCPApplicationFactory")
        .def(py::init<>())
        .def("from_rockit_interface", &fatrop::StageOCPApplicationFactory::from_rockit_interface);

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}