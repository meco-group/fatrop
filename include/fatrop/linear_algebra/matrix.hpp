//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//
#ifndef __fatrop_linear_algebra_matrix_hpp__
#define __fatrop_linear_algebra_matrix_hpp__

/**
 * @file matrix.hpp
 * @brief Provides a C++ interface for working with matrix numerical
 * computations using BLASFEO
 */

#include "fatrop/common/exception.hpp"
#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/blasfeo_wrapper.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>

namespace fatrop
{

    /**
     * @class MatReal
     * @brief Base class for matrix operations.
     *
     * @tparam Derived The derived class implementing specific matrix behavior.
     */
    template <typename Derived> class MatReal
    {
    public:
        /**
         * @brief Accesses the element at the given row and column.
         *
         * @param i Row index of the element.
         * @param j Column index of the element.
         * @return Scalar The value at position (i, j).
         */
        Scalar operator()(const Index i, const Index j) const
        {
            return static_cast<const Derived *>(this)->operator()(i, j);
        }

        /**
         * @brief Gets the number of rows in the matrix.
         *
         * @return Index The number of rows.
         */
        Index m() const { return static_cast<const Derived *>(this)->m(); }

        /**
         * @brief Gets the number of columns in the matrix.
         *
         * @return Index The number of columns.
         */
        Index n() const { return static_cast<const Derived *>(this)->n(); }

        friend MatRealTranspose<Derived> transpose(const MatReal<Derived> &dep)
        {
            return MatRealTranspose<Derived>(dep);
        }

        /**
         * @brief Overloads the << operator for printing the matrix elements.
         *
         * @param os The output stream.
         * @param mat The matrix to print.
         * @return std::ostream& The updated output stream.
         */
        friend std::ostream &operator<<(std::ostream &os, const MatReal<Derived> &mat)
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

    template <typename Dep1> class MatRealTranspose : public MatReal<MatRealTranspose<Dep1>>
    {
    public:
        MatRealTranspose(const MatReal<Dep1> &dep) : dep_(dep) {}
        Scalar operator()(const Index i, const Index j) const { return dep_(j, i); }
        Index m() const { return dep_.n(); }
        Index n() const { return dep_.m(); }

    private:
        const MatReal<Dep1> &dep_;
    };

    /**
     * @class MatRealView1D
     * @brief Represents a 1D view of a matrix (e.g., row, column, or diagonal).
     *
     * @tparam Derived The derived class implementing specific view behavior.
     */
    template <typename Derived> class MatRealView1D : public VecReal<MatRealView1D<Derived>>
    {
    public:
        /**
         * @brief Accesses the element at the given index.
         *
         * @param i Index of the element.
         * @return Scalar& Reference to the value at index i.
         */
        Scalar &operator()(const Index i) { return static_cast<Derived *>(this)->operator()(i); }

        /**
         * @brief Accesses the element at the given index (const version).
         *
         * @param i Index of the element.
         * @return const Scalar The value at index i.
         */
        const Scalar operator()(const Index i) const
        {
            return static_cast<const Derived *>(this)->operator()(i);
        }

        /**
         * @brief Gets the size of the view.
         *
         * @return Index The size of the view.
         */
        Index m() const { return static_cast<const Derived *>(this)->m(); }

        /**
         * @brief Creates a sub-view (block) of the current view.
         *
         * @param size Size of the block.
         * @param start Starting index of the block.
         * @return VecRealBlock<MatRealView1D<Derived>> The sub-view.
         */
        VecRealBlock<MatRealView1D<Derived>> block(const Index size, const Index start) const
        {
            return VecRealBlock<MatRealView1D<Derived>>(*this, size, start);
        }

        /**
         * @brief Assigns values from another vector to this view.
         *
         * @tparam OtherDerived The type of the other vector.
         * @param other The vector to assign from.
         */
        template <typename OtherDerived> void operator=(const VecReal<OtherDerived> &other)
        {
            fatrop_dbg_assert(m() == other.m() && "Vectors must be same size for assignment");
            for (Index i = 0; i < m(); i++)
            {
                (*this)(i) = other(i);
            }
        }

        /**
         * @brief Assigns a scalar value to all elements of the view.
         *
         * @param alpha The scalar value to assign.
         */
        void operator=(const Scalar alpha)
        {
            for (Index i = 0; i < m(); i++)
            {
                (*this)(i) = alpha;
            }
        }
    };

    /**
     * @class MatRealRowView
     * @brief Represents a view of a single row in a matrix.
     */
    class MatRealRowView : public MatRealView1D<MatRealRowView>
    {
    public:
        /**
         * @brief Constructs a MatRealRowView object.
         *
         * @param mat Reference to the matrix.
         * @param row Index of the row to view.
         */
        MatRealRowView(MatRealView &mat, const Index row) : mat_(mat), row_(row) {}

        inline Scalar operator()(const Index i) const;
        inline Scalar &operator()(const Index i);
        inline Index m() const;

        using MatRealView1D<MatRealRowView>::operator=;

    private:
        MatRealView &mat_;
        const Index row_;
    };

    /**
     * @class MatRealColView
     * @brief Represents a view of a single column in a matrix.
     */
    class MatRealColView : public MatRealView1D<MatRealColView>
    {
    public:
        /**
         * @brief Constructs a MatRealColView object.
         *
         * @param mat Reference to the matrix.
         * @param col Index of the column to view.
         */
        MatRealColView(MatRealView &mat, const Index col) : mat_(mat), col_(col) {}

        inline Scalar operator()(const Index i) const;
        inline Scalar &operator()(const Index i);
        inline Index m() const;

        using MatRealView1D<MatRealColView>::operator=;

    private:
        MatRealView &mat_;
        const Index col_;
    };

    /**
     * @class MatRealDiagonalView
     * @brief Represents a view of the diagonal of a matrix.
     */
    class MatRealDiagonalView : public MatRealView1D<MatRealDiagonalView>
    {
    public:
        /**
         * @brief Constructs a MatRealDiagonalView object.
         *
         * @param mat Reference to the matrix.
         */
        MatRealDiagonalView(MatRealView &mat) : mat_(mat) {}

        inline Scalar operator()(const Index i) const;
        inline Scalar &operator()(const Index i);
        inline Index m() const;

        using MatRealView1D<MatRealDiagonalView>::operator=;

    private:
        MatRealView &mat_;
    };

    /**
     * @class MatRealView
     * @brief Represents a numeric matrix with efficient operations.
     */
    class MatRealView : public MatReal<MatRealView>
    {
    public:
        /**
         * @brief Constructs a MatRealView object.
         *
         * @param mat Reference to the allocated matrix.
         * @param m Number of rows.
         * @param n Number of columns.
         * @param ai Row offset.
         * @param aj Column offset.
         */
        MatRealView(MatRealAllocated &mat, const Index m, const Index n, const Index ai,
                    const Index aj)
            : mat_(mat), m_(m), n_(n), ai_(ai), aj_(aj) {};

        /**
         * @brief Accesses the element at the given row and column.
         *
         * @param i Row index of the element.
         * @param j Column index of the element.
         * @return Scalar& Reference to the value at position (i, j).
         */
        inline Scalar &operator()(const Index i, const Index j) const;

        Index m() const { return m_; }
        Index n() const { return n_; }
        Index ai() const { return ai_; }
        Index aj() const { return aj_; }
        inline MAT &mat();
        inline const MAT &mat() const;

        /**
         * @brief Creates a view of a specific row.
         *
         * @param row Index of the row to view.
         * @return MatRealRowView The row view.
         */
        MatRealRowView row(const Index row) { return MatRealRowView(*this, row); }

        /**
         * @brief Creates a view of a specific column.
         *
         * @param col Index of the column to view.
         * @return MatRealColView The column view.
         */
        MatRealColView col(const Index col) { return MatRealColView(*this, col); }

        /**
         * @brief Creates a view of the matrix diagonal.
         *
         * @return MatRealDiagonalView The diagonal view.
         */
        MatRealDiagonalView diagonal() { return MatRealDiagonalView(*this); }

        /**
         * @brief Creates a sub-matrix (block) of the current matrix.
         *
         * @param rows Number of rows in the block.
         * @param cols Number of columns in the block.
         * @param row_start Starting row index of the block.
         * @param col_start Starting column index of the block.
         * @return MatRealView The sub-matrix.
         */
        MatRealView block(const Index rows, const Index cols, const Index row_start,
                          const Index col_start) const
        {
            return MatRealView(mat_, rows, cols, ai_ + row_start, aj_ + col_start);
        }

        /**
         * @brief Assigns a scalar value to all elements of the matrix.
         *
         * @param alpha The scalar value to assign.
         * @return MatRealView& Reference to the modified matrix.
         */
        inline MatRealView &operator=(const Scalar alpha);

        /**
         * @brief Assigns values from another matrix to this matrix.
         *
         * @tparam Derived The type of the other matrix.
         * @param mat_in The matrix to assign from.
         * @return MatRealView& Reference to the modified matrix.
         */
        template <typename Derived> inline MatRealView &operator=(const MatReal<Derived> &mat_in);

        /**
         * @brief Copy assignment operator.
         *
         * @param mat_in The matrix to copy from.
         * @return MatRealView& Reference to the modified matrix.
         */
        inline MatRealView &operator=(const MatRealView &mat_in);

    private:
        MatRealAllocated &mat_;
        const Index m_;
        const Index n_;
        const Index ai_;
        const Index aj_;
    };

    /**
     * @class MatRealAllocated
     * @brief Manages the memory associated with a matrix, providing
     *        allocation and deallocation using BLASFEO functions.
     */
    class MatRealAllocated : public MatRealView
    {
    public:
        /**
         * @brief Constructs a MatRealAllocated object.
         *
         * @param m Number of rows.
         * @param n Number of columns.
         */
        inline MatRealAllocated(const Index m, const Index n);

        /**
         * @brief Deleted copy constructor to prevent unintended copies.
         */
        MatRealAllocated(MatRealAllocated &other) = delete;

        /**
         * @brief Move constructor for MatRealAllocated.
         *
         * @param other The MatRealAllocated object to move from.
         */
        inline MatRealAllocated(MatRealAllocated &&other);

        template <typename Derived>
        MatRealAllocated(const MatReal<Derived> &mat_in) : MatRealAllocated(mat_in.m(), mat_in.n())
        {
            (*this) = mat_in;
        }

        /**
         * @brief Gets a reference to the underlying BLASFEO matrix.
         *
         * @return MAT& Reference to the BLASFEO matrix.
         */
        MAT &mat() { return mat_; }

        /**
         * @brief Gets a const reference to the underlying BLASFEO matrix.
         *
         * @return const MAT& Const reference to the BLASFEO matrix.
         */
        const MAT &mat() const { return mat_; }

        using MatRealView::operator=;

        /**
         * @brief Accesses the element at the given row and column.
         *
         * @param i Row index of the element.
         * @param j Column index of the element.
         * @return Scalar& Reference to the value at position (i, j).
         */
        inline Scalar &operator()(const Index i, const Index j);

        /**
         * @brief Accesses the element at the given row and column (const version).
         *
         * @param i Row index of the element.
         * @param j Column index of the element.
         * @return const Scalar& Const reference to the value at position (i, j).
         */
        inline Scalar operator()(const Index i, const Index j) const;

        MatRealAllocated &operator=(const MatRealAllocated &other)
        {
            static_cast<MatRealView &>(*this) = static_cast<const MatRealView &>(other);
            return *this;
        };

        /**
         * @brief Destructor that frees the allocated matrix memory.
         */
        inline ~MatRealAllocated();

    private:
        MAT mat_;
    };
    // Implementation of MatRealRowView
    Scalar MatRealRowView::operator()(const Index i) const { return mat_(row_, i); }
    Scalar &MatRealRowView::operator()(const Index i) { return mat_(row_, i); }
    Index MatRealRowView::m() const { return mat_.n(); }

    // Implementation of MatrixColumnView
    Scalar MatRealColView::operator()(const Index i) const { return mat_(i, col_); }
    Scalar &MatRealColView::operator()(const Index i) { return mat_(i, col_); }
    Index MatRealColView::m() const { return mat_.m(); }

    // Implementation of MatrixDiagonaView
    Scalar MatRealDiagonalView::operator()(const Index i) const { return mat_(i, i); }
    Scalar &MatRealDiagonalView::operator()(const Index i) { return mat_(i, i); }
    Index MatRealDiagonalView::m() const { return std::min(mat_.m(), mat_.n()); }

    // Implementation of MatRealView
    MAT &MatRealView::mat() { return mat_.mat(); }
    const MAT &MatRealView::mat() const { return mat_.mat(); }
    Scalar &MatRealView::operator()(const Index i, const Index j) const
    {
        return mat_(i + ai_, j + aj_);
    }
    MatRealView &MatRealView::operator=(const Scalar alpha)
    {
        GESE(m(), n(), alpha, &this->mat_.mat(), ai_, aj_);
        return *this;
    }

    template <typename Derived> MatRealView &MatRealView::operator=(const MatReal<Derived> &mat_in)
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
    MatRealView &MatRealView::operator=(const MatRealView &mat_in)
    {
        fatrop_dbg_assert((*this).m() == mat_in.m() && (*this).n() == mat_in.n());
        GECP((*this).m(), (*this).n(), const_cast<MAT *>(&mat_in.mat()), mat_in.ai_, mat_in.aj_,
             &mat_.mat(), ai_, aj_);
        return *this;
    };

    // Implementation of MatRealAllocated methods
    MatRealAllocated::MatRealAllocated(const Index m, const Index n)
        : MatRealView(*this, m, n, 0, 0), mat_()
    {
        ALLOCATE_MAT(m, n, &mat_);
        // zero out
        std::memset(mat_.mem, 0, mat_.memsize * sizeof(char));
    }

    MatRealAllocated::MatRealAllocated(MatRealAllocated &&other)
        : MatRealView(*this, other.m(), other.n(), 0, 0), mat_(other.mat_)
    {
        // Nullify the moved-from object's mat_ to prevent double deletion
        other.mat_.mem = nullptr;
    }

    Scalar &MatRealAllocated::operator()(const Index i, const Index j)
    {
        fatrop_dbg_assert(i >= 0 && i < m() && j >= 0 && j < n());
        return blasfeo_matel_wrap(&mat_, i, j);
    }

    Scalar MatRealAllocated::operator()(const Index i, const Index j) const
    {
        fatrop_dbg_assert(i >= 0 && i < m() && j >= 0 && j < n());
        return blasfeo_matel_wrap(&mat_, i, j);
    }

    MatRealAllocated::~MatRealAllocated()
    {
        if (mat_.mem != nullptr)
            FREE_MAT(&mat());
    }

} // namespace fatrop

#endif // __fatrop_linear_algebra_matrix_hpp__
