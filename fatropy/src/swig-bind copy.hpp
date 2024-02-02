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
        static PySwigObject* get_swig_ptr(PyObject* obj)
        {
            PyObject *this_ptr = PyObject_GetAttrString(obj, "this");
            PySwigObject *swig_obj = reinterpret_cast<PySwigObject *>(this_ptr);
            return swig_obj;
        }
        static T *convert(PyObject *obj_ptr, const std::string &name)
        {
            PySwigObject *swig_obj = get_swig_ptr(obj_ptr);
            std::string py_name = reinterpret_cast<PySwigObject *>(obj_ptr)->ob_base.ob_type->tp_name;
            if (py_name != name)
                throw std::runtime_error("Could not convert of type " + py_name + " to " + name);
            if (swig_obj == nullptr)
                throw std::runtime_error("Could not convert to PySwigObject");
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
                /* Extract PyObject from handle */
                PyObject *source = src.ptr();
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
                pybind11::object ret = attr_(NULL);
                type_cpp* ret_ptr = fatropy::FromPySwig<type_cpp>::convert(ret.ptr(), T::py_name); 
                * ret_ptr = src;
                // PyObject* this_ptr = reinterpret_cast<PyObject*>(fatropy::FromPySwig<type_cpp>::get_swig_ptr(ret.ptr())); 
                Py_IncRef(ret.ptr());
                // Py_IncRef(this_ptr);
                return ret;
            }
        };
    } // namespace PYBIND11_NAMESPACE::detail
}
