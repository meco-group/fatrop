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
#include <iostream>
#include <iomanip>
#include "fatrop/common/exception.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/blasfeo_wrapper.hpp"
#include "fatrop/linear_algebra/vector.hpp"

namespace fatrop
{

    template <typename Derived>
    class Mat
    {
    public:
        Scalar operator()(const Index i, const Index j) const
        {
            return static_cast<const Derived *>(this)->operator()(i, j);
        }

        Index m() const
        {
            return static_cast<const Derived *>(this)->m();
        }

        Index n() const
        {
            return static_cast<const Derived *>(this)->n();
        }

        friend std::ostream &operator<<(std::ostream &os, const Mat<Derived> &mat)
        {
            for (Index i = 0; i < mat.m(); i++)
            {
                for (Index j = 0; j < mat.n(); j++)
                {
                    os << std::setw(12) << std::setprecision(4) << std::fixed << mat(i, j) << " ";
                }
                os << std::endl;
            }
            return os;
        }
    };

    template <typename Derived>
    class MatView1D : public Vec<MatView1D<Derived>>
    {
    public:
        Scalar &operator()(const Index i)
        {
            return static_cast<Derived *>(this)->operator()(i);
        }

        const Scalar operator()(const Index i) const
        {
            return static_cast<const Derived *>(this)->operator()(i);
        }

        Index m() const
        {
            return static_cast<const Derived *>(this)->m();
        }

        VecBlock<MatView1D<Derived>> block(const Index size, const Index start) const
        {
            return VecBlock<MatView1D<Derived>>(*this, size, start);
        }

        template <typename OtherDerived>
        void operator=(const Vec<OtherDerived> &other)
        {
            fatrop_dbg_assert(m() == other.m() && "Vectors must be same size for assignment");
            for (Index i = 0; i < m(); i++)
            {
                (*this)(i) = other(i);
            }
        }

        void operator=(const Scalar alpha)
        {
            for (Index i = 0; i < m(); i++)
            {
                (*this)(i) = alpha;
            }
        }
    };

    class MatRowView : public MatView1D<MatRowView>
    {
    public:
        MatRowView(MatrixNumeric &mat, const Index row) : mat_(mat), row_(row) {}

        Scalar operator()(const Index i) const;
        Scalar &operator()(const Index i);
        Index m() const;

        using MatView1D<MatRowView>::operator=;

    private:
        MatrixNumeric &mat_;
        const Index row_;
    };

    class MatColView : public MatView1D<MatColView>
    {
    public:
        MatColView(MatrixNumeric &mat, const Index col) : mat_(mat), col_(col) {}

        Scalar operator()(const Index i) const;
        Scalar &operator()(const Index i);
        Index m() const;

        using MatView1D<MatColView>::operator=;

    private:
        MatrixNumeric &mat_;
        const Index col_;
    };

    class MatDiagonalView : public MatView1D<MatDiagonalView>
    {
    public:
        MatDiagonalView(MatrixNumeric &mat) : mat_(mat) {}

        Scalar operator()(const Index i) const;
        Scalar &operator()(const Index i);
        Index m() const;

        using MatView1D<MatDiagonalView>::operator=;

    private:
        MatrixNumeric &mat_;
    };

    class MatrixNumeric : public Mat<MatrixNumeric>
    {
    public:
        MatrixNumeric(MatrixAllocated &mat, const Index m, const Index n, const Index ai, const Index aj) : mat_(mat), m_(m), n_(n), ai_(ai), aj_(aj) {};

        Scalar &operator()(const Index i, const Index j) const;

        Index m() const { return m_; }
        Index n() const { return n_; }

        // // View access
        MatRowView row(const Index row) { return MatRowView(*this, row); }
        MatColView col(const Index col) { return MatColView(*this, col); }
        MatDiagonalView diagonal() { return MatDiagonalView(*this); }

        // Block access
        MatrixNumeric block(const Index rows, const Index cols, const Index row_start, const Index col_start) const
        {
            return MatrixNumeric(mat_, rows, cols, ai_ + row_start, aj_ + col_start);
        }

        // Assignment operator
        MatrixNumeric& operator=(const Scalar alpha);
        template <typename Derived>
        MatrixNumeric& operator=(const Mat<Derived> &mat_in);
        MatrixNumeric& operator=(const MatrixNumeric& mat_in){*this = *static_cast<const Mat<MatrixNumeric>*>(&mat_in); return *this;};

    private:
        MatrixAllocated &mat_;
        const Index m_;
        const Index n_;
        const Index ai_;
        const Index aj_;
    };

    /**
     * @class MatrixAllocated
     * @brief Manages the memory associated with providing
     *        allocation and deallocation using BLASFEO functions.
     */
    class MatrixAllocated : public MatrixNumeric
    {
    public:
        MatrixAllocated(const Index m, const Index n);
        MatrixAllocated(MatrixAllocated &other) = delete;
        MatrixAllocated(MatrixAllocated &&other);

        MAT &mat() { return mat_; }
        const MAT &mat() const { return mat_; }

        using MatrixNumeric::operator=;

        // Element access
        Scalar &operator()(const Index i, const Index j);
        const Scalar &operator()(const Index i, const Index j) const;

        ~MatrixAllocated();

    private:
        MAT mat_;
    };
    // Implementation of MatRowView
    Scalar MatRowView::operator()(const Index i) const { return mat_(row_, i); }
    Scalar &MatRowView::operator()(const Index i) { return mat_(row_, i); }
    Index MatRowView::m() const { return mat_.n(); }

    // Implementation of MatrixColumnView
    Scalar MatColView::operator()(const Index i) const { return mat_(i, col_); }
    Scalar &MatColView::operator()(const Index i) { return mat_(i, col_); }
    Index MatColView::m() const { return mat_.m(); }

    // Implementation of MatrixDiagonaView
    Scalar MatDiagonalView::operator()(const Index i) const { return mat_(i, i); }
    Scalar &MatDiagonalView::operator()(const Index i) { return mat_(i, i); }
    Index MatDiagonalView::m() const { return std::min(mat_.m(), mat_.n()); }

    // Implementation of MatrixNumeric
    Scalar &MatrixNumeric::operator()(const Index i, const Index j) const
    {
        return mat_(i + ai_, j + aj_);
    }
    MatrixNumeric& MatrixNumeric::operator=(const Scalar alpha)
    {
        for (Index i = 0; i < m_; i++)
        {
            for (Index j = 0; j < n_; j++)
            {
                (*this)(i, j) = alpha;
            }
        }
        return *this;
    }

    template <typename Derived>
    MatrixNumeric& MatrixNumeric::operator=(const Mat<Derived> &mat_in)
    {
        fatrop_dbg_assert((*this).m() == mat_in.m());
        fatrop_dbg_assert((*this).n() == mat_in.n());
        for (Index i = 0; i < m_; i++)
        {
            for (Index j = 0; j < n_; j++)
            {
                (*this)(i, j) = mat_in(i, j);
            }
        }
        return *this;
    }

    // Implementation of MatrixAllocated methods
    MatrixAllocated::MatrixAllocated(const Index m, const Index n) : MatrixNumeric(*this, m, n, 0, 0), mat_()
    {
        ALLOCATE_MAT(m, n, &mat_);
    }

    MatrixAllocated::MatrixAllocated(MatrixAllocated &&other)
        : MatrixNumeric(*this, other.m(), other.n(), 0, 0), mat_(other.mat_)
    {
        // Nullify the moved-from object's mat_ to prevent double deletion
        other.mat_.mem = nullptr;
    }

    Scalar &MatrixAllocated::operator()(const Index i, const Index j)
    {
        fatrop_dbg_assert(i >= 0 && i < m() && j >= 0 && j < n());
        return MATEL(&mat_, i, j);
    }

    const Scalar &MatrixAllocated::operator()(const Index i, const Index j) const
    {
        fatrop_dbg_assert(i >= 0 && i < m() && j >= 0 && j < n());
        return MATEL(&mat_, i, j);
    }

    MatrixAllocated::~MatrixAllocated()
    {
        FREE_MAT(&mat());
    }

} // namespace fatrop

#endif // __fatrop_linear_algebra_matrix_hpp__
