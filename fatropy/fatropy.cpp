#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "ocp/StageOCPApplication.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
namespace fatropy
{
    struct MatBind : public std::vector<double>
    {
        using std::vector<double>::vector;
        MatBind(const py::array_t<double> &x)
        {
            auto buf = x.request();
            auto ptr = static_cast<double *>(buf.ptr);
            this->resize(buf.size);
            if (buf.ndim == 1)
            {
                std::copy(ptr, ptr + buf.size, this->begin());
            }
            else
            {
                const int m = buf.shape[0];
                const int n = buf.shape[1];
                const int stride_0 = buf.strides[0] / buf.itemsize;
                const int stride_1 = buf.strides[1] / buf.itemsize;
                // copy to fortran format
                for (int j = 0; j < n; j++)
                    for (int i = 0; i < m; i++)
                        (*this)[i + j * m] = ptr[stride_0 * i + stride_1 * j];
            }
        };
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
        .def("set_initial", &fatrop::NLPApplication::set_initial)
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

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}