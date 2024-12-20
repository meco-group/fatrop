// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_linear_algebra_matrix_hpp__
#define __fatrop_linear_algebra_matrix_hpp__

/**
 * @file matrix.hpp
 * @brief Provides a C++ interface for working with matrix numerical
 * computations using BLASFEO (Basic Linear Algebra Subprograms For Embedded
 * Optimization).
 */

#include <cmath>
#include "fatrop/common/exception.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/blasfeo_wrapper.hpp"

namespace fatrop
{

    /**
     * @class MatrixAllocated
     * @brief Manages the memory associated with providing
     *        allocation and deallocation using BLASFEO functions.
     */
    class MatrixAllocated
    {
    public:
        MatrixAllocated(const Index m, const Index n);
        MatrixAllocated(MatrixAllocated &other) = delete;
        MatrixAllocated(MatrixAllocated &&other);


        MAT &mat() { return mat_; }
        const MAT &mat() const { return mat_; }

        Index m() const { return m_; }
        Index n() const { return n_; }

        // Element access
        Scalar &operator()(const Index i, const Index j);
        const Scalar &operator()(const Index i, const Index j) const;

        ~MatrixAllocated();

    private:
        MAT mat_;
        const Index m_;
        const Index n_;
    };

    // Implementation of MatrixAllocated methods
    MatrixAllocated::MatrixAllocated(const Index m, const Index n):m_(m), n_(n)
    {
        ALLOCATE_MAT(m, n, &mat_);
    }

    MatrixAllocated::MatrixAllocated(MatrixAllocated &&other)
        : mat_(other.mat_), m_(other.m_), n_(other.n_)
    {
        // Nullify the moved-from object's mat_ to prevent double deletion
        other.mat_.mem = nullptr;
    }

    Scalar &MatrixAllocated::operator()(const Index i, const Index j)
    {
        fatrop_dbg_assert(i >= 0 && i < m_ && j >= 0 && j < n_);
        return MATEL(&mat_, i, j);
    }

    const Scalar &MatrixAllocated::operator()(const Index i, const Index j) const
    {
        fatrop_dbg_assert(i >= 0 && i < m_ && j >= 0 && j < n_);
        return MATEL(&mat_, i, j);
    }

    MatrixAllocated::~MatrixAllocated()
    {
        FREE_MAT(&mat());
    }

} // namespace fatrop

#endif // __fatrop_linear_algebra_matrix_hpp__
