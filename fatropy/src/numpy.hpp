#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>
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