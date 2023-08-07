#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <casadi/casadi.hpp>
#include <vector>
namespace fatropy
{
    namespace py = pybind11;
    struct PySwigObject
    {
        PyObject_HEAD
        void *ptr;
        const char* decr;
    };
    template <typename T>
    struct FromPySwig
    {
        static T *convert(const py::object &obj)
        {
            // todo: figure out how this actually works, it is inspired by eigenpy
            PyObject *obj_ptr = obj.ptr();
            PyObject* this_ptr = PyObject_GetAttrString(obj_ptr, "this");
            PySwigObject *swig_obj = reinterpret_cast<PySwigObject *>(this_ptr);
            if(this_ptr == nullptr)
                throw std::runtime_error("Could not convert to PySwigObject");
            auto ret = reinterpret_cast<T *>(swig_obj->ptr);
            return ret;
        }
    };
    
    void print_function_nametest(casadi::Function &f)
    {
        std::cout << "Function pointer: ";
        std::cout << &f << std::endl;
        std::cout << "Function name: ";
        std::cout << f.name() << std::endl;
        std::cout << "number of inputs = " << f.n_in() << std::endl;
    };
}