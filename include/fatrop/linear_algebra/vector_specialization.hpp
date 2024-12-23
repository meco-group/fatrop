//
// Copyright (c) Lander Vanroye, KU Leuven
//

#ifndef __fatrop_linear_algebra_vector_specialization_hpp__
#define __fatrop_linear_algebra_vector_specialization_hpp__

#include "fatrop/common/exception.hpp"
/**
 * @file vector_specializations.hpp
 * @brief implements specializations for redirecting vector assignment operations to blasfeo
 * kernels. These specializations are not striclty necessary as their functionality is already
 * included in "vector.hpp", but they might improve performance.
 */

namespace fatrop
{
    /**
     * @class VecRealViewPlusVecRealView
     * @brief Represents the addition of two `VecRealView` vectors.
     */
    class VecRealViewPlusVecRealView : public VecOperationSpecialization<VecRealViewPlusVecRealView>
    {
    public:
        VecRealViewPlusVecRealView(VecRealView &a, VecRealView &b) : a_(a), b_(b)
        {
            fatrop_dbg_assert(a.m() == b.m() && "Vector sizes must match");
        };
        Scalar operator()(const Index i) const { return a_(i) + b_(i); }
        Index m() const { return a_.m(); }
        friend inline VecRealViewPlusVecRealView operator+(VecRealView &a, VecRealView &b);
        VecRealView &a() const { return a_; }
        VecRealView &b() const { return b_; }

    private:
        VecRealView &a_;
        VecRealView &b_;
    };

    /**
     * @class VecRealViewTimesScalar
     * @brief Represents the scaling of a `VecRealView` vector by a scalar.
     */
    class VecRealViewTimesScalar : public VecOperationSpecialization<VecRealViewTimesScalar>
    {
    public:
        VecRealViewTimesScalar(VecRealView &a, const Scalar alpha) : a_(a), alpha_(alpha) {};
        Scalar operator()(const Index i) const { return alpha_ * a_(i); }
        Index m() const { return a_.m(); }
        friend inline VecRealViewTimesScalar operator*(const Scalar alpha, VecRealView &a);
        friend inline VecRealViewTimesScalar operator*(VecRealView &a, const Scalar alpha);
        friend inline VecRealViewTimesScalar operator-(VecRealView &a);
        VecRealView &a() const { return a_; }
        Scalar alpha() const { return alpha_; }

    private:
        VecRealView &a_;
        const Scalar alpha_;
    };

    /**
     * @class VecRealViewTimesScalarPlusVecRealViewTimesScalar
     * @brief Represents a linear combination of two scaled `VecRealView` vectors.
     */

    class VecRealViewTimesScalarPlusVecRealViewTimesScalar
        : public VecOperationSpecialization<VecRealViewTimesScalarPlusVecRealViewTimesScalar>
    {
    public:
        VecRealViewTimesScalarPlusVecRealViewTimesScalar(VecRealView &a, const Scalar alpha,
                                                       VecRealView &b, const Scalar beta)
            : a_(a), alpha_(alpha), b_(b), beta_(beta)
        {
            fatrop_dbg_assert(a.m() == b.m() && "Vector sizes must match");
        };
        Scalar operator()(const Index i) const { return alpha_ * a_(i) + beta_ * b_(i); }
        Index m() const { return a_.m(); }
        friend inline VecRealViewTimesScalarPlusVecRealViewTimesScalar
        operator+(const VecRealViewTimesScalar &a, const VecRealViewTimesScalar &b);
        VecRealView &a() const { return a_; }
        Scalar alpha() const { return alpha_; }
        VecRealView &b() const { return b_; }
        Scalar beta() const { return beta_; }

    private:
        VecRealView &a_;
        const Scalar alpha_;
        VecRealView &b_;
        const Scalar beta_;
    };

    /**
     * @class VecRealViewPlusVecRealViewTimesScalar
     * @brief Represents the addition of a `VecRealView` vector and a scaled
     * `VecRealView` vector.
     */
    class VecRealViewPlusVecRealViewTimesScalar
        : public VecOperationSpecialization<VecRealViewPlusVecRealViewTimesScalar>
    {
    public:
        VecRealViewPlusVecRealViewTimesScalar(VecRealView &a, VecRealView &b, const Scalar alpha)
            : a_(a), b_(b), alpha_(alpha)
        {
            fatrop_dbg_assert(a.m() == b.m() && "Vector sizes must match");
        };
        Scalar operator()(const Index i) const { return a_(i) + alpha_ * b_(i); }
        Index m() const { return a_.m(); }
        friend inline VecRealViewPlusVecRealViewTimesScalar operator+(const VecRealViewTimesScalar &a,
                                                             VecRealView &b);
        friend inline VecRealViewPlusVecRealViewTimesScalar operator+(VecRealView &a,
                                                             const VecRealViewTimesScalar &b);
        VecRealView &a() const { return a_; }
        VecRealView &b() const { return b_; }
        Scalar alpha() const { return alpha_; }

    private:
        VecRealView &a_;
        VecRealView &b_;
        const Scalar alpha_;
    };

    /**
     * @class VecRealViewTimesVecRealView
     * @brief Represents the element-wise multiplication of two `VecRealView`
     * vectors.
     */
    class VecRealViewTimesVecRealView : public VecOperationSpecialization<VecRealViewTimesVecRealView>
    {
    public:
        VecRealViewTimesVecRealView(VecRealView &a, VecRealView &b) : a_(a), b_(b)
        {
            fatrop_dbg_assert(a.m() == b.m() && "Vector sizes must match");
        };
        Scalar operator()(const Index i) const { return a_(i) * b_(i); }
        Index m() const { return a_.m(); }
        friend inline VecRealViewTimesVecRealView operator*(VecRealView &a, VecRealView &b);
        VecRealView &a() const { return a_; }
        VecRealView &b() const { return b_; }

    private:
        VecRealView &a_;
        VecRealView &b_;
    };

    // operator overloading for VecRealView specializations - blasfeo kernels
    template <>
    inline VecRealView &VecRealView::operator=(
        VecOperationSpecialization<VecRealViewPlusVecRealView> &&vecnumericplusvecnumeric)
    {
        auto vv = vecnumericplusvecnumeric.derived();
        fatrop_dbg_assert(this->m() == vv.a().m() && "Vector sizes must match");
        VecRealView &a = vv.a();
        VecRealView &b = vv.b();
        AXPY(m_, 1.0, &a.vec(), a.ai(), &b.vec(), b.ai(), &this->vec(), this->ai());
        return *this;
    }

    template <>
    inline VecRealView &
    VecRealView::operator=(VecOperationSpecialization<VecRealViewTimesScalar> &&vecnumerictimesscalar)
    {
        auto vv = vecnumerictimesscalar.derived();
        fatrop_dbg_assert(this->m() == vv.a().m() && "Vector sizes must match");
        VecRealView &a = vv.a();
        Scalar alpha = vv.alpha();
        VECCPSC(m_, alpha, &a.vec(), a.ai(), &this->vec(), this->ai());
        return *this;
    }

    template <>
    inline VecRealView &
    VecRealView::operator=(VecOperationSpecialization<VecRealViewPlusVecRealViewTimesScalar>
                              &&vecnumericplusvecnumerictimesscalar)
    {
        auto vv = vecnumericplusvecnumerictimesscalar.derived();
        fatrop_dbg_assert(this->m() == vv.a().m() && "Vector sizes must match");
        VecRealView &a = vv.a();
        VecRealView &b = vv.b();
        Scalar alpha = vv.alpha();
        AXPY(m_, alpha, &b.vec(), b.ai(), &a.vec(), a.ai(), &this->vec(), this->ai());
        return *this;
    }

    template <>
    inline VecRealView &
    VecRealView::operator=(VecOperationSpecialization<VecRealViewTimesScalarPlusVecRealViewTimesScalar>
                              &&vecnumerictimesscalarplusvecnumerictimesscalar)
    {
        auto vv = vecnumerictimesscalarplusvecnumerictimesscalar.derived();
        fatrop_dbg_assert(this->m() == vv.a().m() && "Vector sizes must match");
        VecRealView &a = vv.a();
        Scalar alpha = vv.alpha();
        VecRealView &b = vv.b();
        Scalar beta = vv.beta();
        AXPBY(m_, alpha, &a.vec(), a.ai(), beta, &b.vec(), b.ai(), &this->vec(), this->ai());
        return *this;
    }

    template <>
    inline VecRealView &VecRealView::operator=(
        VecOperationSpecialization<VecRealViewTimesVecRealView> &&vecnumerictimesscalar)
    {
        auto vv = vecnumerictimesscalar.derived();
        fatrop_dbg_assert(this->m() == vv.a().m() && "Vector sizes must match");
        VecRealView &a = vv.a();
        VecRealView &b = vv.b();
        VECMUL(m_, &b.vec(), b.ai(), &a.vec(), a.ai(), &this->vec(), this->ai());
        return *this;
    }

    /**
     * @brief Addition operator for VecRealView and VecRealViewTimesScalar.
     */
    VecRealViewPlusVecRealViewTimesScalar operator+(VecRealView &a, const VecRealViewTimesScalar &b)
    {
        return VecRealViewPlusVecRealViewTimesScalar(a, b.a(), b.alpha());
    }

    /**
     * @brief Addition operator for VecRealViewTimesScalar and VecRealView.
     */
    VecRealViewPlusVecRealViewTimesScalar operator+(const VecRealViewTimesScalar &a, VecRealView &b)
    {
        return b + a;
    }

    /**
     * @brief Addition operator for two VecRealViewTimesScalar objects.
     */
    VecRealViewTimesScalarPlusVecRealViewTimesScalar operator+(const VecRealViewTimesScalar &a,
                                                             const VecRealViewTimesScalar &b)
    {
        return VecRealViewTimesScalarPlusVecRealViewTimesScalar(a.a(), a.alpha(), b.a(), b.alpha());
    }

    /**
     * @brief Element-wise multiplication operator for two VecRealView objects.
     */
    VecRealViewTimesVecRealView operator*(VecRealView &a, VecRealView &b)
    {
        return VecRealViewTimesVecRealView(a, b);
    }

    /**
     * @brief Addition operator for two VecRealView objects.
     */
    VecRealViewPlusVecRealView operator+(VecRealView &a, VecRealView &b)
    {
        return VecRealViewPlusVecRealView(a, b);
    }

    /**
     * @brief Scalar multiplication operator for VecRealView.
     */
    VecRealViewTimesScalar operator*(const Scalar alpha, VecRealView &a)
    {
        return VecRealViewTimesScalar(a, alpha);
    }

    /**
     * @brief Scalar multiplication operator for VecRealView.
     */
    VecRealViewTimesScalar operator*(VecRealView &a, const Scalar alpha) { return alpha * a; }

    /**
     * @brief Unary minus operator for VecRealView.
     */
    VecRealViewTimesScalar operator-(VecRealView &a) { return -1. * a; }

} // namespace fatrop

#endif // __fatrop_vector_specialization_hpp__
