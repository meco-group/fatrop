#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <casadi/casadi.hpp>
#include <vector>
namespace py = pybind11;
namespace fatropy
{
    struct PySwigObject
    {
        PyObject_HEAD void *ptr;
        const char *decr;
    };
    template <typename T>
    struct FromPySwig
    {
        static T *convert(const py::object &obj, const std::string &name)
        {
            // todo: figure out how this actually works, it is inspired by eigenpy
            PyObject *obj_ptr = obj.ptr();
            return convert(obj_ptr, name);
        }
        static T *convert(PyObject *obj_ptr, const std::string &name)
        {
            PyObject *this_ptr = PyObject_GetAttrString(obj_ptr, "this");
            PySwigObject *swig_obj = reinterpret_cast<PySwigObject *>(this_ptr);
            std::string py_name = reinterpret_cast<PySwigObject *>(obj_ptr)->ob_base.ob_type->tp_name;
            if (py_name != name)
                throw std::runtime_error("Could not convert of type " + py_name + " to " + name);
            if (this_ptr == nullptr)
                throw std::runtime_error("Could not convert to PySwigObject");
            Py_DECREF(this_ptr);
            auto ret = reinterpret_cast<T *>(swig_obj->ptr);
            return ret;
        }
    };
}
namespace PYBIND11_NAMESPACE
{
    namespace detail
    {
        template <typename T>
        struct swig_type_caster
        {
        public:
            using type_cpp = typename T::type;
            // using name_cpp = T::name;
            PYBIND11_TYPE_CASTER(type_cpp, const_name(T::py_name));

            bool load(handle src, bool)
            {
                // check if handle is a list 
                if(py::isinstance<py::list>(src)) return false;
                pybind11::module_ cspy_ = pybind11::module_::import(T::module);
                pybind11::object attr_ = cspy_.attr(T::py_name);
                pybind11::object src2(attr_(src));
                /* Extract PyObject from handle */
                PyObject *source = src2.ptr();
                /* Try converting into a Python integer value */
                PyObject *tmp = source;
                if (!tmp)
                    return false;
                /* Now try to convert into a C++ int */
                // value.long_value = PyLong_AsLong(tmp);
                value = *fatropy::FromPySwig<type_cpp>::convert(source, T::py_name);
                // Py_DECREF(tmp);
                /* Ensure return code was OK (to avoid out-of-range errors etc) */
                return !PyErr_Occurred();
                // return false;
            }

            static handle cast(type_cpp src, return_value_policy /* policy */, handle /* parent */)
            {
                pybind11::module_ cspy_ = pybind11::module_::import(T::module);
                pybind11::object attr_ = cspy_.attr(T::py_name);
                pybind11::object ret(attr_());
                *fatropy::FromPySwig<type_cpp>::convert(ret, T::py_name) = src;
                return ret.release();
            }
        };
    } // namespace PYBIND11_NAMESPACE::detail
}
