//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//
#ifndef __fatrop_linear_algebra_vector_hpp__
#define __fatrop_linear_algebra_vector_hpp__

/**
 * @file vector.hpp
 * @brief Provides a C++ interface for working with vectorized numerical
 * computations using BLASFEO
 */

#include "fatrop/common/exception.hpp"
#include "fatrop/context/context.hpp"
#include "fatrop/linear_algebra/blasfeo_wrapper.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include <cmath>
#include <cstring>
#include <iomanip>
#include <vector>

namespace fatrop
{
    template <typename Derived> class VecReal
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
         * @return VecRealBlock<Derived> The sub-block of the vector.
         */
        VecRealBlock<Derived> block(Index size, Index start) const
        {
            return VecRealBlock<Derived>(*static_cast<const Derived *>(this), size, start);
        }

        const Derived &derived() const { return *static_cast<const Derived *>(this); }

        bool is_zero() const
        {
            for (Index i = 0; i < m(); i++)
            {
                if ((*this)(i) != 0.)
                {
                    return false;
                }
            }
            return true;
        }

        // Various mathematical operations defined as friend functions
        // They allow mathematical operations like sum, norms, and transformations on
        // vectors.
        friend Scalar sum(const VecReal<Derived> &vec)
        {
            Scalar ret = 0;
            for (Index i = 0; i < vec.m(); i++)
            {
                ret += vec(i);
            }
            return ret;
        }

        friend Scalar norm_inf(const VecReal<Derived> &vec)
        {
            Scalar ret = 0;
            for (Index i = 0; i < vec.m(); i++)
            {
                ret = std::max(ret, std::abs(vec(i)));
            }
            return ret;
        }

        friend Scalar norm_l1(const VecReal<Derived> &vec)
        {
            Scalar ret = 0;
            for (Index i = 0; i < vec.m(); i++)
            {
                ret += std::abs(vec(i));
            }
            return ret;
        }

        friend Scalar norm_l2(const VecReal<Derived> &vec)
        {
            Scalar ret = 0;
            for (Index i = 0; i < vec.m(); i++)
            {
                ret += vec(i) * vec(i);
            }
            return std::sqrt(ret);
        }

        template <typename Dep2>
        friend Scalar dot(const VecReal<Derived> &a, const VecReal<Dep2> &b)
        {
            Scalar ret = 0;
            for (Index i = 0; i < a.m(); i++)
            {
                ret += a(i) * b(i);
            }
            return ret;
        }

        friend VecRealAbs<Derived> abs(const VecReal<Derived> &a) { return VecRealAbs<Derived>(a); }

        friend VecRealLog<Derived> log(const VecReal<Derived> &a) { return VecRealLog<Derived>(a); }

        friend VecRealExp<Derived> exp(const VecReal<Derived> &a) { return VecRealExp<Derived>(a); }

        friend VecRealSin<Derived> sin(const VecReal<Derived> &a) { return VecRealSin<Derived>(a); }

        friend VecRealCos<Derived> cos(const VecReal<Derived> &a) { return VecRealCos<Derived>(a); }

        template <typename Dep2>
        friend VecRealMin<Derived, Dep2> min(const VecReal<Derived> &a, const VecReal<Dep2> &b)
        {
            return VecRealMin<Derived, Dep2>(a, b);
        }

        template <typename Dep2>
        friend VecRealMax<Derived, Dep2> max(const VecReal<Derived> &a, const VecReal<Dep2> &b)
        {
            return VecRealMax<Derived, Dep2>(a, b);
        }

        template <typename IfElseOp, typename Dep2>
        friend VecRealIfElse<IfElseOp, Derived, Dep2>
        if_else(const IfElseOp &if_else_op, const VecReal<Derived> &a, const VecReal<Dep2> &b)
        {
            return VecRealIfElse<IfElseOp, Derived, Dep2>(if_else_op, a, b);
        }
        // Operator overloading functions of Expressions
        template <typename Dep2>
        friend VecRealPlusVecReal<Derived, Dep2> operator+(const VecReal<Derived> &a,
                                                           const VecReal<Dep2> &b)
        {
            return VecRealPlusVecReal<Derived, Dep2>(a, b);
        }

        template <typename Dep2>
        friend VecRealMinusVecReal<Derived, Dep2> operator-(const VecReal<Derived> &a,
                                                            const VecReal<Dep2> &b)
        {
            return VecRealMinusVecReal<Derived, Dep2>(a, b);
        }

        friend VecRealTimesScalar<Derived> operator*(const Scalar alpha, const VecReal<Derived> &a)
        {
            return VecRealTimesScalar<Derived>(a, alpha);
        }

        friend VecRealTimesScalar<Derived> operator*(const VecReal<Derived> &a, const Scalar alpha)
        {
            return alpha * a;
        }

        friend VecRealTimesScalar<Derived> operator-(VecReal<Derived> &a) { return -1. * a; }

        template <typename Dep2>
        friend VecRealTimesVecReal<Derived, Dep2> operator*(const VecReal<Derived> &a,
                                                            const VecReal<Dep2> &b)
        {
            return VecRealTimesVecReal<Derived, Dep2>(a, b);
        }

        template <typename Dep2>
        friend VecRealDivVecReal<Derived, Dep2> operator/(const VecReal<Derived> &a,
                                                          const VecReal<Dep2> &b)
        {
            return VecRealDivVecReal<Derived, Dep2>(a, b);
        }

        friend ScalarDivVecReal<Derived> operator/(const Scalar alpha, const VecReal<Derived> &a)
        {
            return ScalarDivVecReal<Derived>(alpha, a);
        }

        /**
         * @brief Overloads the << operator for printing the vector elements.
         *
         * @param os The output stream.
         * @param vec The vector to print.
         * @return std::ostream& The updated output stream.
         */
        friend std::ostream &operator<<(std::ostream &os, const VecReal<Derived> &vec)
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
     * @class VecRealBlock
     * @brief Represents a sub-segment (block) of a vector.
     *
     * @tparam Dep1 Type of the parent vector.
     */
    template <typename Dep1> class VecRealBlock : public VecReal<VecRealBlock<Dep1>>
    {
    public:
        VecRealBlock(const VecReal<Dep1> &a, const Index m, const Index ai)
            : a(a.derived()), m_(m), ai_(ai) {};
        Scalar operator()(const Index i) const { return a(ai_ + i); }
        Index m() const { return m_; }

    private:
        const Dep1 a;
        const Index m_;
        const Index ai_;
    };

    // Other vector operations are similarly defined below.
    /**
     * @brief Represents an element-wise conditional vector operation.
     */
    template <typename IfElseOp, typename Dep1, typename Dep2>
    class VecRealIfElse : public VecReal<VecRealIfElse<IfElseOp, Dep1, Dep2>>
    {
    public:
        VecRealIfElse(const IfElseOp &if_else_op, const VecReal<Dep1> &a, const VecReal<Dep2> &b)
            : if_else_op(if_else_op), a(a.derived()), b(b.derived())
        {
        }
        Scalar operator()(const Index i) const { return if_else_op(i) ? a(i) : b(i); }
        Index m() const { return a.m(); }

    private:
        const Dep1 a;
        const Dep2 b;
        const IfElseOp &if_else_op;
    };

    /**
     * @brief Represents element-wise minimum of two vectors.
     */
    template <typename Dep1, typename Dep2>
    class VecRealMin : public VecReal<VecRealMin<Dep1, Dep2>>
    {
    public:
        VecRealMin(const VecReal<Dep1> &a, const VecReal<Dep2> &b) : a(a.derived()), b(b.derived())
        {
            fatrop_dbg_assert(a.m() == b.m() && "Vector sizes must match for min operation");
        }
        Scalar operator()(const Index i) const { return std::min(a(i), b(i)); }
        Index m() const { return a.m(); }

    private:
        const Dep1 a;
        const Dep2 b;
    };

    /**
     * @brief Represents element-wise maximum of two vectors.
     */
    template <typename Dep1, typename Dep2>
    class VecRealMax : public VecReal<VecRealMax<Dep1, Dep2>>
    {
    public:
        VecRealMax(const VecReal<Dep1> &a, const VecReal<Dep2> &b) : a(a.derived()), b(b.derived())
        {
            fatrop_dbg_assert(a.m() == b.m() && "Vector sizes must match for max operation");
        }
        Scalar operator()(const Index i) const { return std::max(a(i), b(i)); }
        Index m() const { return a.m(); }

    private:
        const Dep1 a;
        const Dep2 b;
    };

    // specialization for std::vector<bool>
    template <typename Dep1, typename Dep2>
    class VecRealIfElse<std::vector<bool>, Dep1, Dep2>
        : public VecReal<VecRealIfElse<std::vector<bool>, Dep1, Dep2>>
    {
    public:
        VecRealIfElse(const std::vector<bool> &if_else_op, const VecReal<Dep1> &a,
                      const VecReal<Dep2> &b)
            : if_else_op(if_else_op), a(a.derived()), b(b.derived())
        {
        }
        Scalar operator()(const Index i) const { return if_else_op[i] ? a(i) : b(i); }
        Index m() const { return a.m(); }

    private:
        const Dep1 a;
        const Dep2 b;
        const std::vector<bool> &if_else_op;
    };

    /**
     * @brief Represents a vector filled with a single scalar value.
     */
    class VecRealScalar : public VecReal<VecRealScalar>
    {
    public:
        explicit VecRealScalar(const Index m, const Scalar alpha) : m_(m), alpha(alpha) {};
        Scalar operator()(const Index i) const { return alpha; }
        Index m() const { return m_; }

    private:
        const Index m_;
        const Scalar alpha;
    };

    /**
     * @brief Represents element-wise addition of two vectors.
     */
    template <typename Dep1, typename Dep2>
    class VecRealPlusVecReal : public VecReal<VecRealPlusVecReal<Dep1, Dep2>>
    {
    public:
        VecRealPlusVecReal(const VecReal<Dep1> &a, const VecReal<Dep2> &b) : a(a.derived()), b(b.derived())
        {
            fatrop_dbg_assert(a.m() == b.m() && "Vector sizes must match for addition");
        };
        Scalar operator()(const Index i) const { return a(i) + b(i); }
        Index m() const { return a.m(); }

    private:
        const Dep1 a;
        const Dep2 b;
    };

    /**
     * @brief Represents element-wise subtraction of two vectors.
     */
    template <typename Dep1, typename Dep2>
    class VecRealMinusVecReal : public VecReal<VecRealMinusVecReal<Dep1, Dep2>>
    {
    public:
        VecRealMinusVecReal(const VecReal<Dep1> &a, const VecReal<Dep2> &b) : a(a.derived()), b(b.derived())
        {
            fatrop_dbg_assert(a.m() == b.m() && "Vector sizes must match for subtraction");
        };
        Scalar operator()(const Index i) const { return a(i) - b(i); }
        Index m() const { return a.m(); }

    private:
        const Dep1 a;
        const Dep2 b;
    };

    /**
     * @brief Represents scalar multiplication of a vector.
     */
    template <typename Dep1> class VecRealTimesScalar : public VecReal<VecRealTimesScalar<Dep1>>
    {
    public:
        VecRealTimesScalar(const VecReal<Dep1> &a, const Scalar alpha) : a(a.derived()), alpha(alpha) {};
        Scalar operator()(const Index i) const { return alpha * a(i); }
        Index m() const { return a.m(); }

    private:
        const Dep1 a;
        const Scalar alpha;
    };

    /**
     * @brief Represents element-wise multiplication of two vectors.
     */
    template <typename Dep1, typename Dep2>
    class VecRealTimesVecReal : public VecReal<VecRealTimesVecReal<Dep1, Dep2>>
    {
    public:
        VecRealTimesVecReal(const VecReal<Dep1> &a, const VecReal<Dep2> &b) : a(a.derived()), b(b.derived())
        {
            fatrop_dbg_assert(a.m() == b.m() &&
                              "Vector sizes must match for element-wise multiplication");
        };
        Scalar operator()(const Index i) const { return a(i) * b(i); }
        Index m() const { return a.m(); }

    private:
        const Dep1 a;
        const Dep2 b;
    };

    /**
     * @brief Represents element-wise division of two vectors.
     */
    template <typename Dep1, typename Dep2>
    class VecRealDivVecReal : public VecReal<VecRealDivVecReal<Dep1, Dep2>>
    {
    public:
        VecRealDivVecReal(const VecReal<Dep1> &a, const VecReal<Dep2> &b) : a(a.derived()), b(b.derived())
        {
            fatrop_dbg_assert(a.m() == b.m() &&
                              "Vector sizes must match for element-wise division");
        };
        Scalar operator()(const Index i) const { return a(i) / b(i); }
        Index m() const { return a.m(); }

    private:
        const Dep1 a;
        const Dep2 b;
    };

    /**
     * @brief Represents scalar division by a vector.
     */
    template <typename Dep1> class ScalarDivVecReal : public VecReal<ScalarDivVecReal<Dep1>>
    {
    public:
        ScalarDivVecReal(const Scalar alpha, const VecReal<Dep1> &a) : alpha(alpha), a(a.derived()) {};
        Scalar operator()(const Index i) const { return alpha / a(i); }
        Index m() const { return a.m(); }

    private:
        const Scalar alpha;
        const Dep1 a;
    };

    /**
     * @brief Represents element-wise absolute value of a vector.
     */
    template <typename Dep1> class VecRealAbs : public VecReal<VecRealAbs<Dep1>>
    {
    public:
        VecRealAbs(const VecReal<Dep1> &a) : a(a.derived()) {};
        Scalar operator()(const Index i) const { return std::abs(a(i)); }
        Index m() const { return a.m(); }

    private:
        const Dep1 a;
    };

    /**
     * @brief Represents element-wise natural logarithm of a vector.
     */
    template <typename Dep1> class VecRealLog : public VecReal<VecRealLog<Dep1>>
    {
    public:
        VecRealLog(const VecReal<Dep1> &a) : a(a.derived()) {};
        Scalar operator()(const Index i) const { return std::log(a(i)); }
        Index m() const { return a.m(); }

    private:
        const Dep1 a;
    };

    /**
     * @brief Represents element-wise exponential of a vector.
     */
    template <typename Dep1> class VecRealExp : public VecReal<VecRealExp<Dep1>>
    {
    public:
        VecRealExp(const VecReal<Dep1> &a) : a(a.derived()) {};
        Scalar operator()(const Index i) const { return std::exp(a(i)); }
        Index m() const { return a.m(); }

    private:
        const Dep1 a;
    };

    /**
     * @brief Represents element-wise sine of a vector.
     */
    template <typename Dep1> class VecRealSin : public VecReal<VecRealSin<Dep1>>
    {
    public:
        VecRealSin(const VecReal<Dep1> &a) : a(a.derived()) {};
        Scalar operator()(const Index i) const { return std::sin(a(i)); }
        Index m() const { return a.m(); }

    private:
        const Dep1 a;
    };

    /**
     * @brief Represents element-wise cosine of a vector.
     */
    template <typename Dep1> class VecRealCos : public VecReal<VecRealCos<Dep1>>
    {
    public:
        VecRealCos(const VecReal<Dep1> &a) : a(a.derived()) {};
        Scalar operator()(const Index i) const { return std::cos(a(i)); }
        Index m() const { return a.m(); }

    private:
        const Dep1 a;
    };

    /**
     * @brief Base class for vector operation specializations.
     */
    template <typename Derived> class VecOperationSpecialization : public VecReal<Derived>
    {
    public:
        const Derived &derived() const { return static_cast<const Derived &>(*this); }
    };

    /**
     * @class VecRealView
     * @brief Represents a numerically stored vector with efficient operations
     *        backed by the BLASFEO library.
     */
    class VecRealView : public VecReal<VecRealView>
    {
    public:
        inline VecRealView(VecRealAllocated &vec, const Index m, const Index ai);
        inline VecRealView(VecRealAllocated &vec);
        VecRealView(const VecReal<VecRealView> &vec) : VecRealView(vec.derived()) {};
        // blasfeo_dvecse
        inline VecRealView &operator=(const Scalar alpha);
        // Assignment operators to handle various vector operations efficiently
        // by redirecting to blasfeo kernels
        // assignment from VecOperationSpecialization
        template <typename Derived>
        inline VecRealView &operator=(VecOperationSpecialization<Derived> &&vec_in);
        // assignment from general VecNumeric expression
        template <typename Derived> inline VecRealView &operator=(const VecReal<Derived> &vec_in);
        VecRealView &operator=(const VecRealView &vec_in)
        {
            *this = *static_cast<const VecReal<VecRealView> *>(&vec_in);
            return *this;
        };
        /**
         * @brief Accessor for elements of the vector.
         */
        inline Scalar &operator()(const Index i) const;
        /**
         * @brief Creates a sub-vector (block) of the current vector.
         */
        VecRealView block(Index size, Index start) const
        {
            return VecRealView(vec_, size, ai_ + start);
        }
        inline VEC &vec();
        inline const VEC &vec() const;
        inline Index m() const;
        inline Index ai() const;
        Scalar *data() { return m() > 0 ? &this->operator()(0) : nullptr; }
        const Scalar *data() const { return m() > 0 ? &this->operator()(0) : nullptr; }

    private:
        VecRealAllocated &vec_;
        const Index m_;
        const Index ai_;
    };

    /**
     * @class VecRealAllocated
     * @brief Manages the memory associated with a VecRealView, providing
     *        allocation and deallocation using BLASFEO functions.
     */
    class VecRealAllocated : public VecRealView
    {
    public:
        /**
         * @brief Constructs a VecRealAllocated object with the given size.
         */
        VecRealAllocated(const Index m) : VecRealView(*this, m, 0), m_(m)
        {
            ALLOCATE_VEC(m, &vec_);
            std::memset(vec_.mem, 0, vec_.memsize * sizeof(char));
        }
        template <typename Derived>
        VecRealAllocated(const VecReal<Derived> &vec) : VecRealAllocated(vec.m())
        {
            *this = vec;
        }
        VecRealAllocated(const VecRealAllocated &other) = delete;
        /**
         * @brief Move constructor for VecRealAllocated.
         */
        VecRealAllocated(VecRealAllocated &&other)
            : VecRealView(*this, other.m(), 0), vec_(other.vec()), m_(other.m())
        {
            // Nullify the moved-from object's vec_ to prevent double deletion
            other.vec_.mem = nullptr;
        }
        using VecRealView::operator=;
        VEC &vec() { return vec_; }
        const VEC &vec() const { return vec_; }

        Index m() const { return m_; }
        /**
         * @brief Accessor for mutable elements of the vector.
         */
        Scalar &operator()(const Index i)
        {
            fatrop_dbg_assert(i >= 0 && i < m_);
            return blasfeo_vecel_wrap(&vec(), i);
        }
        /**
         * @brief Accessor for const elements of the vector.
         */
        Scalar operator()(const Index i) const
        {
            fatrop_dbg_assert(i >= 0 && i < m_);
            VEC veccp = vec();
            return blasfeo_vecel_wrap(&vec(), i);
        }

        VecRealAllocated &operator=(const VecRealAllocated &other)
        {
            static_cast<VecRealView &>(*this) = static_cast<const VecRealView &>(other);
            return *this;
        };
        /**
         * @brief Destructor that frees the allocated vector memory.
         */
        ~VecRealAllocated()
        {
            if (vec_.mem != nullptr)
                FREE_VEC(&vec());
        }

    private:
        VEC vec_;
        const Index m_;
    };

    // implementation
    VecRealView::VecRealView(VecRealAllocated &vec, const Index m, const Index ai)
        : vec_(vec), m_(m), ai_(ai) {};
    VecRealView ::VecRealView(VecRealAllocated &vec) : vec_(vec), m_(vec.m()), ai_(0) {};
    template <typename Derived>
    inline VecRealView &VecRealView::operator=(const VecReal<Derived> &vec_in)
    {
        fatrop_dbg_assert(m() == vec_in.m() && "Vectors must be same size for asignment");
        for (Index i = 0; i < m(); i++)
        {
            this->operator()(i) = vec_in(i);
        }
        return *this;
    }
    VecRealView &VecRealView::operator=(const Scalar alpha)
    {
        VECSE(m(), alpha, &vec(), ai());
        return *this;
    }

    Scalar &VecRealView::operator()(const Index i) const
    {
        fatrop_dbg_assert(i >= 0 && i < m_);
        return vec_(ai() + i);
    }
    VEC &VecRealView::vec() { return vec_.vec(); }
    const VEC &VecRealView::vec() const { return vec_.vec(); }
    Index VecRealView::m() const { return m_; }
    Index VecRealView::ai() const { return ai_; }

} // namespace fatrop

// add specializations
#include "vector_specialization.hpp"

#endif // __fatrop_linear_algebra_vector_hpp__
