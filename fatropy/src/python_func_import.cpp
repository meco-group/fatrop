#include "python_func_import.hpp"
#include "expose-spectool.hpp"
#include "pybind11/embed.h"
namespace fatropy
{
    casadi::Function get_func_from_py(const std::string &file_name, const std::string &func_name)
    {
        py::scoped_interpreter guard{};
        // pybind11 import func_name from file_name (absolute path)
        py::module m = py::module_::import(file_name.c_str());
        py::object func = m.attr(func_name.c_str());
        // call func
        py::object ret = func();
        return ret.cast<casadi::Function>();
    }
    void donothgin(){};
}
