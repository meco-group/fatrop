//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_linear_algebra_vector_hpp__
#define __fatrop_linear_algebra_vector_hpp__

/**
 * @file vector.hpp
 * @brief Provides a C++ interface for working with vectorized numerical
 * computations using BLASFEO (Basic Linear Algebra Subprograms For Embedded
 * Optimization).
 */

#include <cmath>
#include <iomanip>
#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/blasfeo_wrapper.hpp"
#include "fatrop/common/exception.hpp"
#include "fatrop/linear_algebra/fwd.hpp"

namespace fatrop
{
    template <typename Derived> class Vec
    {
    public:
        /**
         * @brief Accesses the element at the given index.
         *
         * @param i Index of the element.
         * @return Scalar The value at index i.
         */
        Scalar operator()(const Index i) const
        {
            return static_cast<const Derived *>(this)->operator()(i);
        }

        /**
         * @brief Gets the size (number of elements) of the vector.
         *
         * @return int The size of the vector.
         */
        Index m() const { return static_cast<const Derived *>(this)->m(); }

        /**
         * @brief Extracts a sub-block (segment) of the vector.
         *
         * @param size Size of the block.
         * @param start Starting index of the block.
         * @return VecBlock<Derived> The sub-block of the vector.
         */
        VecBlock<Derived> block(Index size, Index start) const
        {
            return VecBlock<Derived>(*static_cast<const Derived *>(this), size, start);
        }

        // Various mathematical operations defined as friend functions
        // They allow mathematical operations like sum, norms, and transformations on
        // vectors.
        friend Scalar sum(const Vec<Derived> &vec)
        {
            Scalar ret = 0;
            for (Index i = 0; i < vec.m(); i++)
            {
                ret += vec(i);
            }
            return ret;
        }

        friend Scalar norm_inf(const Vec<Derived> &vec)
        {
            Scalar ret = 0;
            for (Index i = 0; i < vec.m(); i++)
            {
                ret = std::max(ret, std::abs(vec(i)));
            }
            return ret;
        }

        friend Scalar norm_l1(const Vec<Derived> &vec)
        {
            Scalar ret = 0;
            for (Index i = 0; i < vec.m(); i++)
            {
                ret += std::abs(vec(i));
            }
            return ret;
        }

        friend Scalar norm_l2(const Vec<Derived> &vec)
        {
            Scalar ret = 0;
            for (Index i = 0; i < vec.m(); i++)
            {
                ret += vec(i) * vec(i);
            }
            return std::sqrt(ret);
        }

        template <typename Dep2> friend Scalar dot(const Vec<Derived> &a, const Vec<Dep2> &b)
        {
            Scalar ret = 0;
            for (Index i = 0; i < a.m(); i++)
            {
                ret += a(i) * b(i);
            }
            return ret;
        }

        friend VecAbs<Derived> abs(const Vec<Derived> &a) { return VecAbs<Derived>(a); }

        friend VecLog<Derived> log(const Vec<Derived> &a) { return VecLog<Derived>(a); }

        friend VecExp<Derived> exp(const Vec<Derived> &a) { return VecExp<Derived>(a); }

        friend VecSin<Derived> sin(const Vec<Derived> &a) { return VecSin<Derived>(a); }

        friend VecCos<Derived> cos(const Vec<Derived> &a) { return VecCos<Derived>(a); }

        template <typename IfElseOp, typename Dep2>
        friend VecIfElse<IfElseOp, Derived, Dep2> if_else(const IfElseOp &if_else_op,
                                                          const Vec<Derived> &a, const Vec<Dep2> &b)
        {
            return VecIfElse<IfElseOp, Derived, Dep2>(if_else_op, a, b);
        }
        // Operator overloading functions of Expressions
        template <typename Dep2>
        friend VecPlusVec<Derived, Dep2> operator+(const Vec<Derived> &a, const Vec<Dep2> &b)
        {
            return VecPlusVec<Derived, Dep2>(a, b);
        }

        template <typename Dep2>
        friend VecMinusVec<Derived, Dep2> operator-(const Vec<Derived> &a, const Vec<Dep2> &b)
        {
            return VecMinusVec<Derived, Dep2>(a, b);
        }

        friend VecTimesScalar<Derived> operator*(const Scalar alpha, const Vec<Derived> &a)
        {
            return VecTimesScalar<Derived>(a, alpha);
        }

        friend VecTimesScalar<Derived> operator*(const Vec<Derived> &a, const Scalar alpha)
        {
            return alpha * a;
        }

        friend VecTimesScalar<Derived> operator-(Vec<Derived> &a) { return -1. * a; }

        template <typename Dep2>
        friend VecTimesVec<Derived, Dep2> operator*(const Vec<Derived> &a, const Vec<Dep2> &b)
        {
            return VecTimesVec<Derived, Dep2>(a, b);
        }

        template <typename Dep2>
        friend VecDivVec<Derived, Dep2> operator/(const Vec<Derived> &a, const Vec<Dep2> &b)
        {
            return VecDivVec<Derived, Dep2>(a, b);
        }

        friend ScalarDivVec<Derived> operator/(const Scalar alpha, const Vec<Derived> &a)
        {
            return ScalarDivVec<Derived>(alpha, a);
        }

        /**
         * @brief Overloads the << operator for printing the vector elements.
         *
         * @param os The output stream.
         * @param vec The vector to print.
         * @return std::ostream& The updated output stream.
         */
        friend std::ostream &operator<<(std::ostream &os, const Vec<Derived> &vec)
        {
            for (Index i = 0; i < vec.m(); i++)
            {
                os << std::setw(12) << std::setprecision(4) << vec(i) << " ";
            }
            return os; // Return the stream to allow chaining (e.g., cout << obj1 <<
                       // obj2)
        }
    };

    /**
     * @class VecBlock
     * @brief Represents a sub-segment (block) of a vector.
     *
     * @tparam Dep1 Type of the parent vector.
     */
    template <typename Dep1> class VecBlock : public Vec<VecBlock<Dep1>>
    {
    public:
        VecBlock(const Vec<Dep1> &a, const Index m, const Index ai) : a(a), m_(m), ai_(ai) {};
        Scalar operator()(const Index i) const { return a(ai_ + i); }
        Index m() const { return m_; }

    private:
        const Vec<Dep1> &a;
        const Index m_;
        const Index ai_;
    };

    // Other vector operations are similarly defined below.
    /**
     * @brief Represents an element-wise conditional vector operation.
     */
    template <typename IfElseOp, typename Dep1, typename Dep2>
    class VecIfElse : public Vec<VecIfElse<IfElseOp, Dep1, Dep2>>
    {
    public:
        VecIfElse(const IfElseOp &if_else_op, const Vec<Dep1> &a, const Vec<Dep2> &b)
            : if_else_op(if_else_op), a(a), b(b)
        {
        }
        Scalar operator()(const Index i) const { return if_else_op(i) ? a(i) : b(i); }
        Index m() const { return a.m(); }

    private:
        const Vec<Dep1> &a;
        const Vec<Dep2> &b;
        const IfElseOp &if_else_op;
    };

    /**
     * @brief Represents a vector filled with a single scalar value.
     */
    class VecScalar : public Vec<VecScalar>
    {
    public:
        explicit VecScalar(const Index m, const Scalar alpha) : m_(m), alpha(alpha) {};
        Scalar operator()(const Index i) const { return alpha; }
        Index m() const { return m_; }

    private:
        const Index m_;
        const Scalar alpha;
    };

    /**
     * @brief Represents element-wise addition of two vectors.
     */
    template <typename Dep1, typename Dep2> class VecPlusVec : public Vec<VecPlusVec<Dep1, Dep2>>
    {
    public:
        VecPlusVec(const Vec<Dep1> &a, const Vec<Dep2> &b) : a(a), b(b)
        {
            fatrop_dbg_assert(a.m() == b.m() && "Vector sizes must match for addition");
        };
        Scalar operator()(const Index i) const { return a(i) + b(i); }
        Index m() const { return a.m(); }

    private:
        const Vec<Dep1> &a;
        const Vec<Dep2> &b;
    };

    /**
     * @brief Represents element-wise subtraction of two vectors.
     */
    template <typename Dep1, typename Dep2> class VecMinusVec : public Vec<VecMinusVec<Dep1, Dep2>>
    {
    public:
        VecMinusVec(const Vec<Dep1> &a, const Vec<Dep2> &b) : a(a), b(b)
        {
            fatrop_dbg_assert(a.m() == b.m() && "Vector sizes must match for subtraction");
        };
        Scalar operator()(const Index i) const { return a(i) - b(i); }
        Index m() const { return a.m(); }

    private:
        const Vec<Dep1> &a;
        const Vec<Dep2> &b;
    };

    /**
     * @brief Represents scalar multiplication of a vector.
     */
    template <typename Dep1> class VecTimesScalar : public Vec<VecTimesScalar<Dep1>>
    {
    public:
        VecTimesScalar(const Vec<Dep1> &a, const Scalar alpha) : a(a), alpha(alpha) {};
        Scalar operator()(const Index i) const { return alpha * a(i); }
        Index m() const { return a.m(); }

    private:
        const Vec<Dep1> &a;
        const Scalar alpha;
    };

    /**
     * @brief Represents element-wise multiplication of two vectors.
     */
    template <typename Dep1, typename Dep2> class VecTimesVec : public Vec<VecTimesVec<Dep1, Dep2>>
    {
    public:
        VecTimesVec(const Vec<Dep1> &a, const Vec<Dep2> &b) : a(a), b(b)
        {
            fatrop_dbg_assert(a.m() == b.m() &&
                              "Vector sizes must match for element-wise multiplication");
        };
        Scalar operator()(const Index i) const { return a(i) * b(i); }
        Index m() const { return a.m(); }

    private:
        const Vec<Dep1> &a;
        const Vec<Dep2> &b;
    };

    /**
     * @brief Represents element-wise division of two vectors.
     */
    template <typename Dep1, typename Dep2> class VecDivVec : public Vec<VecDivVec<Dep1, Dep2>>
    {
    public:
        VecDivVec(const Vec<Dep1> &a, const Vec<Dep2> &b) : a(a), b(b)
        {
            fatrop_dbg_assert(a.m() == b.m() &&
                              "Vector sizes must match for element-wise division");
        };
        Scalar operator()(const Index i) const { return a(i) / b(i); }
        Index m() const { return a.m(); }

    private:
        const Vec<Dep1> &a;
        const Vec<Dep2> &b;
    };

    /**
     * @brief Represents scalar division by a vector.
     */
    template <typename Dep1> class ScalarDivVec : public Vec<ScalarDivVec<Dep1>>
    {
    public:
        ScalarDivVec(const Scalar alpha, const Vec<Dep1> &a) : alpha(alpha), a(a) {};
        Scalar operator()(const Index i) const { return alpha / a(i); }
        Index m() const { return a.m(); }

    private:
        const Scalar alpha;
        const Vec<Dep1> &a;
    };

    /**
     * @brief Represents element-wise absolute value of a vector.
     */
    template <typename Dep1> class VecAbs : public Vec<VecAbs<Dep1>>
    {
    public:
        VecAbs(const Vec<Dep1> &a) : a(a) {};
        Scalar operator()(const Index i) const { return std::abs(a(i)); }
        Index m() const { return a.m(); }

    private:
        const Vec<Dep1> &a;
    };

    /**
     * @brief Represents element-wise natural logarithm of a vector.
     */
    template <typename Dep1> class VecLog : public Vec<VecLog<Dep1>>
    {
    public:
        VecLog(const Vec<Dep1> &a) : a(a) {};
        Scalar operator()(const Index i) const { return std::log(a(i)); }
        Index m() const { return a.m(); }

    private:
        const Vec<Dep1> &a;
    };

    /**
     * @brief Represents element-wise exponential of a vector.
     */
    template <typename Dep1> class VecExp : public Vec<VecExp<Dep1>>
    {
    public:
        VecExp(const Vec<Dep1> &a) : a(a) {};
        Scalar operator()(const Index i) const { return std::exp(a(i)); }
        Index m() const { return a.m(); }

    private:
        const Vec<Dep1> &a;
    };

    /**
     * @brief Represents element-wise sine of a vector.
     */
    template <typename Dep1> class VecSin : public Vec<VecSin<Dep1>>
    {
    public:
        VecSin(const Vec<Dep1> &a) : a(a) {};
        Scalar operator()(const Index i) const { return std::sin(a(i)); }
        Index m() const { return a.m(); }

    private:
        const Vec<Dep1> &a;
    };

    /**
     * @brief Represents element-wise cosine of a vector.
     */
    template <typename Dep1> class VecCos : public Vec<VecCos<Dep1>>
    {
    public:
        VecCos(const Vec<Dep1> &a) : a(a) {};
        Scalar operator()(const Index i) const { return std::cos(a(i)); }
        Index m() const { return a.m(); }

    private:
        const Vec<Dep1> &a;
    };

    /**
     * @brief Base class for vector operation specializations.
     */
    template <typename Derived> class VecOperationSpecialization : public Vec<Derived>
    {
    public:
        const Derived &derived() const { return static_cast<const Derived &>(*this); }
    };

    /**
     * @class VecNumeric
     * @brief Represents a numerically stored vector with efficient operations
     *        backed by the BLASFEO library.
     */
    class VecNumeric : public Vec<VecNumeric>
    {
    public:
        inline VecNumeric(VecAllocated &vec, const Index m, const Index ai);
        inline VecNumeric(VecAllocated &vec);
        // blasfeo_dvecse
        inline VecNumeric &operator=(const Scalar alpha);
        // Assignment operators to handle various vector operations efficiently
        // by redirecting to blasfeo kernels
        // assignment from VecOperationSpecialization
        template <typename Derived>
        inline VecNumeric &operator=(VecOperationSpecialization<Derived> &&vec_in);
        // assignment from general Vec expression
        template <typename Derived> inline VecNumeric &operator=(const Vec<Derived> &vec_in);
        VecNumeric &operator=(const VecNumeric &vec_in)
        {
            *this = *static_cast<const Vec<VecNumeric> *>(&vec_in);
            return *this;
        };
        /**
         * @brief Accessor for elements of the vector.
         */
        inline Scalar &operator()(const Index i) const;
        /**
         * @brief Creates a sub-vector (block) of the current vector.
         */
        VecNumeric block(Index size, Index start) const
        {
            return VecNumeric(vec_, size, ai_ + start);
        }
        inline VEC &vec();
        inline Index m() const;
        inline Index ai() const;

    private:
        VecAllocated &vec_;
        const Index m_;
        const Index ai_;
    };

    /**
     * @class VecAllocated
     * @brief Manages the memory associated with a VecNumeric, providing
     *        allocation and deallocation using BLASFEO functions.
     */
    class VecAllocated : public VecNumeric
    {
    public:
        /**
         * @brief Constructs a VecAllocated object with the given size.
         */
        VecAllocated(const Index m) : VecNumeric(*this, m, 0), m_(m) { ALLOCATE_VEC(m, &vec_); }
        VecAllocated(VecAllocated & /*other*/) = delete;
        /**
         * @brief Move constructor for VecAllocated.
         */
        VecAllocated(VecAllocated &&other)
            : VecNumeric(*this, other.m(), 0), vec_(other.vec()), m_(other.m())
        {
            // Nullify the moved-from object's vec_ to prevent double deletion
            other.vec_.mem = nullptr;
        }
        using VecNumeric::operator=;
        VEC &vec() { return vec_; }
        const VEC &vec() const { return vec_; }

        Index m() const { return m_; }
        /**
         * @brief Accessor for mutable elements of the vector.
         */
        Scalar &operator()(const Index i)
        {
            fatrop_dbg_assert(i >= 0 && i < m_);
            return VECEL(&vec(), i);
        }
        /**
         * @brief Accessor for const elements of the vector.
         */
        const Scalar &operator()(const Index i) const
        {
            fatrop_dbg_assert(i >= 0 && i < m_);
            VEC veccp = vec();
            return VECEL(&vec(), i);
        }
        /**
         * @brief Destructor that frees the allocated vector memory.
         */
        ~VecAllocated() { FREE_VEC(&vec()); }

    private:
        VEC vec_;
        const Index m_;
    };

    // implementation
    VecNumeric::VecNumeric(VecAllocated &vec, const Index m, const Index ai)
        : vec_(vec), m_(m), ai_(ai) {};
    VecNumeric ::VecNumeric(VecAllocated &vec) : vec_(vec), m_(vec.m()), ai_(0) {};
    template <typename Derived> inline VecNumeric &VecNumeric::operator=(const Vec<Derived> &vec_in)
    {
        fatrop_dbg_assert(m() == vec_in.m() && "Vectors must be same size for asignment");
        for (Index i = 0; i < m(); i++)
        {
            this->operator()(i) = vec_in(i);
        }
        return *this;
    }
    VecNumeric &VecNumeric::operator=(const Scalar alpha)
    {
        VECSE(m(), alpha, &vec(), ai());
        return *this;
    }

    Scalar &VecNumeric::operator()(const Index i) const
    {
        fatrop_dbg_assert(i >= 0 && i < m_);
        return vec_(ai() + i);
    }
    VEC &VecNumeric::vec() { return vec_.vec(); }
    Index VecNumeric::m() const { return m_; }
    Index VecNumeric::ai() const { return ai_; }

} // namespace fatrop

// add specializations
#include "vector-specialization.hpp"

#endif // __fatrop_linear_algebra_vector_hpp__
