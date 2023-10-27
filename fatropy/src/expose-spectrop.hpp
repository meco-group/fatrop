#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "ocp/StageOCPApplication.hpp"
#include "numpy.hpp"
#include "swig-bind.hpp"
#include "spectrop.hpp"

namespace py = pybind11;

struct ExposeSpectrop
{
    static void expose(py::module &m)
    {
        py::class_<fatrop::spectrop::Ocp>(m, "Ocp").def(py::init<>());
    }
};