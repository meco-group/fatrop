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
#include <iomanip>
#include <iostream>

namespace fatrop
{

    /**
     * @class Mat
     * @brief Base class for matrix operations.
     *
     * @tparam Derived The derived class implementing specific matrix behavior.
     */
    template <typename Derived> class Mat
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

        /**
         * @brief Overloads the << operator for printing the matrix elements.
         *
         * @param os The output stream.
         * @param mat The matrix to print.
         * @return std::ostream& The updated output stream.
         */
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

    /**
     * @class MatView1D
     * @brief Represents a 1D view of a matrix (e.g., row, column, or diagonal).
     *
     * @tparam Derived The derived class implementing specific view behavior.
     */
    template <typename Derived> class MatView1D : public Vec<MatView1D<Derived>>
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
         * @return VecBlock<MatView1D<Derived>> The sub-view.
         */
        VecBlock<MatView1D<Derived>> block(const Index size, const Index start) const
        {
            return VecBlock<MatView1D<Derived>>(*this, size, start);
        }

        /**
         * @brief Assigns values from another vector to this view.
         *
         * @tparam OtherDerived The type of the other vector.
         * @param other The vector to assign from.
         */
        template <typename OtherDerived> void operator=(const Vec<OtherDerived> &other)
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
     * @class MatRowView
     * @brief Represents a view of a single row in a matrix.
     */
    class MatRowView : public MatView1D<MatRowView>
    {
    public:
        /**
         * @brief Constructs a MatRowView object.
         *
         * @param mat Reference to the matrix.
         * @param row Index of the row to view.
         */
        MatRowView(MatrixNumeric &mat, const Index row) : mat_(mat), row_(row) {}

        inline Scalar operator()(const Index i) const;
        inline Scalar &operator()(const Index i);
        inline Index m() const;

        using MatView1D<MatRowView>::operator=;

    private:
        MatrixNumeric &mat_;
        const Index row_;
    };

    /**
     * @class MatColView
     * @brief Represents a view of a single column in a matrix.
     */
    class MatColView : public MatView1D<MatColView>
    {
    public:
        /**
         * @brief Constructs a MatColView object.
         *
         * @param mat Reference to the matrix.
         * @param col Index of the column to view.
         */
        MatColView(MatrixNumeric &mat, const Index col) : mat_(mat), col_(col) {}

        inline Scalar operator()(const Index i) const;
        inline Scalar &operator()(const Index i);
        inline Index m() const;

        using MatView1D<MatColView>::operator=;

    private:
        MatrixNumeric &mat_;
        const Index col_;
    };

    /**
     * @class MatDiagonalView
     * @brief Represents a view of the diagonal of a matrix.
     */
    class MatDiagonalView : public MatView1D<MatDiagonalView>
    {
    public:
        /**
         * @brief Constructs a MatDiagonalView object.
         *
         * @param mat Reference to the matrix.
         */
        MatDiagonalView(MatrixNumeric &mat) : mat_(mat) {}

        inline Scalar operator()(const Index i) const;
        inline Scalar &operator()(const Index i);
        inline Index m() const;

        using MatView1D<MatDiagonalView>::operator=;

    private:
        MatrixNumeric &mat_;
    };

    /**
     * @class MatrixNumeric
     * @brief Represents a numeric matrix with efficient operations.
     */
    class MatrixNumeric : public Mat<MatrixNumeric>
    {
    public:
        /**
         * @brief Constructs a MatrixNumeric object.
         *
         * @param mat Reference to the allocated matrix.
         * @param m Number of rows.
         * @param n Number of columns.
         * @param ai Row offset.
         * @param aj Column offset.
         */
        MatrixNumeric(MatrixAllocated &mat, const Index m, const Index n, const Index ai,
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

        /**
         * @brief Creates a view of a specific row.
         *
         * @param row Index of the row to view.
         * @return MatRowView The row view.
         */
        MatRowView row(const Index row) { return MatRowView(*this, row); }

        /**
         * @brief Creates a view of a specific column.
         *
         * @param col Index of the column to view.
         * @return MatColView The column view.
         */
        MatColView col(const Index col) { return MatColView(*this, col); }

        /**
         * @brief Creates a view of the matrix diagonal.
         *
         * @return MatDiagonalView The diagonal view.
         */
        MatDiagonalView diagonal() { return MatDiagonalView(*this); }

        /**
         * @brief Creates a sub-matrix (block) of the current matrix.
         *
         * @param rows Number of rows in the block.
         * @param cols Number of columns in the block.
         * @param row_start Starting row index of the block.
         * @param col_start Starting column index of the block.
         * @return MatrixNumeric The sub-matrix.
         */
        MatrixNumeric block(const Index rows, const Index cols, const Index row_start,
                            const Index col_start) const
        {
            return MatrixNumeric(mat_, rows, cols, ai_ + row_start, aj_ + col_start);
        }

        /**
         * @brief Assigns a scalar value to all elements of the matrix.
         *
         * @param alpha The scalar value to assign.
         * @return MatrixNumeric& Reference to the modified matrix.
         */
        inline MatrixNumeric &operator=(const Scalar alpha);

        /**
         * @brief Assigns values from another matrix to this matrix.
         *
         * @tparam Derived The type of the other matrix.
         * @param mat_in The matrix to assign from.
         * @return MatrixNumeric& Reference to the modified matrix.
         */
        template <typename Derived> inline MatrixNumeric &operator=(const Mat<Derived> &mat_in);

        /**
         * @brief Copy assignment operator.
         *
         * @param mat_in The matrix to copy from.
         * @return MatrixNumeric& Reference to the modified matrix.
         */
        inline MatrixNumeric &operator=(const MatrixNumeric &mat_in);

    private:
        MatrixAllocated &mat_;
        const Index m_;
        const Index n_;
        const Index ai_;
        const Index aj_;
    };

    /**
     * @class MatrixAllocated
     * @brief Manages the memory associated with a matrix, providing
     *        allocation and deallocation using BLASFEO functions.
     */
    class MatrixAllocated : public MatrixNumeric
    {
    public:
        /**
         * @brief Constructs a MatrixAllocated object.
         *
         * @param m Number of rows.
         * @param n Number of columns.
         */
        inline MatrixAllocated(const Index m, const Index n);

        /**
         * @brief Deleted copy constructor to prevent unintended copies.
         */
        MatrixAllocated(MatrixAllocated &other) = delete;

        /**
         * @brief Move constructor for MatrixAllocated.
         *
         * @param other The MatrixAllocated object to move from.
         */
        inline MatrixAllocated(MatrixAllocated &&other);

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

        using MatrixNumeric::operator=;

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
        inline const Scalar &operator()(const Index i, const Index j) const;

        /**
         * @brief Destructor that frees the allocated matrix memory.
         */
        inline ~MatrixAllocated();

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
    MatrixNumeric &MatrixNumeric::operator=(const Scalar alpha)
    {
        GESE(m(), n(), alpha, &this->mat_.mat(), ai_, aj_);
        return *this;
    }

    template <typename Derived> MatrixNumeric &MatrixNumeric::operator=(const Mat<Derived> &mat_in)
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
    MatrixNumeric &MatrixNumeric::operator=(const MatrixNumeric &mat_in)
    {
        fatrop_dbg_assert((*this).m() == mat_in.m() && (*this).n() == mat_in.n());
        GECP((*this).m(), (*this).n(), &mat_.mat(), ai_, aj_, &mat_in.mat_.mat(), mat_in.ai_,
             mat_in.aj_);
        return *this;
    };

    // Implementation of MatrixAllocated methods
    MatrixAllocated::MatrixAllocated(const Index m, const Index n)
        : MatrixNumeric(*this, m, n, 0, 0), mat_()
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

    MatrixAllocated::~MatrixAllocated() { FREE_MAT(&mat()); }

} // namespace fatrop

#endif // __fatrop_linear_algebra_matrix_hpp__
