//
// Copyright (c) Lander Vanroye, KU Leuven
//

#ifndef __fatrop_linear_algebra_vector_specialization_hpp__
#define __fatrop_linear_algebra_vector_specialization_hpp__

#include "fatrop/common/exception.hpp"

namespace fatrop
{
    /**
     * @class VecNumericPlusVecNumeric
     * @brief Represents the addition of two `VecNumeric` vectors.
     */
    class VecNumericPlusVecNumeric
        : public VecOperationSpecialization<VecNumericPlusVecNumeric>
    {
    public:
        VecNumericPlusVecNumeric(VecNumeric &a, VecNumeric &b) : a_(a), b_(b) {
            fatrop_dbg_assert(a.m() == b.m() && "Vector sizes must match");
        };
        Scalar operator()(const Index i) const { return a_(i) + b_(i); }
        Index m() const { return a_.m(); }
        friend VecNumericPlusVecNumeric operator+(VecNumeric &a, VecNumeric &b);
        VecNumeric& a() const { return a_; }
        VecNumeric& b() const { return b_; }
    private:
        VecNumeric &a_;
        VecNumeric &b_;
    };

    /**
     * @class VecNumericTimesScalar
     * @brief Represents the scaling of a `VecNumeric` vector by a scalar.
     */
    class VecNumericTimesScalar
        : public VecOperationSpecialization<VecNumericTimesScalar>
    {
    public:
        VecNumericTimesScalar(VecNumeric &a, const Scalar alpha)
            : a_(a), alpha_(alpha) {};
        Scalar operator()(const Index i) const { return alpha_ * a_(i); }
        Index m() const { return a_.m(); }
        friend VecNumericTimesScalar operator*(const Scalar alpha, VecNumeric &a);
        friend VecNumericTimesScalar operator*(VecNumeric &a, const Scalar alpha);
        friend VecNumericTimesScalar operator-(VecNumeric &a);
        VecNumeric& a() const { return a_; }
        Scalar alpha() const { return alpha_; }
    private:
        VecNumeric &a_;
        const Scalar alpha_;
    };

    /**
     * @class VecNumericTimesScalarPlusVecNumericTimesScalar
     * @brief Represents a linear combination of two scaled `VecNumeric` vectors.
     */

    class VecNumericTimesScalarPlusVecNumericTimesScalar
        : public VecOperationSpecialization<
              VecNumericTimesScalarPlusVecNumericTimesScalar>
    {
    public:
        VecNumericTimesScalarPlusVecNumericTimesScalar(VecNumeric &a,
                                                       const Scalar alpha,
                                                       VecNumeric &b,
                                                       const Scalar beta)
            : a_(a), alpha_(alpha), b_(b), beta_(beta) {
            fatrop_dbg_assert(a.m() == b.m() && "Vector sizes must match");
        };
        Scalar operator()(const Index i) const { return alpha_ * a_(i) + beta_ * b_(i); }
        Index m() const { return a_.m(); }
        friend VecNumericTimesScalarPlusVecNumericTimesScalar operator+(
            const VecNumericTimesScalar &a, const VecNumericTimesScalar &b);
        VecNumeric& a() const { return a_; }
        Scalar alpha() const { return alpha_; }
        VecNumeric& b() const { return b_; }
        Scalar beta() const { return beta_; }
    private:
        VecNumeric &a_;
        const Scalar alpha_;
        VecNumeric &b_;
        const Scalar beta_;
    };

    /**
     * @class VecNumericPlusVecNumericTimesScalar
     * @brief Represents the addition of a `VecNumeric` vector and a scaled
     * `VecNumeric` vector.
     */
    class VecNumericPlusVecNumericTimesScalar
        : public VecOperationSpecialization<VecNumericPlusVecNumericTimesScalar>
    {
    public:
        VecNumericPlusVecNumericTimesScalar(VecNumeric &a, VecNumeric &b,
                                            const Scalar alpha)
            : a_(a), b_(b), alpha_(alpha) {
            fatrop_dbg_assert(a.m() == b.m() && "Vector sizes must match");
        };
        Scalar operator()(const Index i) const { return a_(i) + alpha_ * b_(i); }
        Index m() const { return a_.m(); }
        friend VecNumericPlusVecNumericTimesScalar operator+(
            const VecNumericTimesScalar &a, VecNumeric &b);
        friend VecNumericPlusVecNumericTimesScalar operator+(
            VecNumeric &a, const VecNumericTimesScalar &b);
        VecNumeric& a() const { return a_; }
        VecNumeric& b() const { return b_; }
        Scalar alpha() const { return alpha_; }
    private:
        VecNumeric &a_;
        VecNumeric &b_;
        const Scalar alpha_;
    };

    /**
     * @class VecNumericTimesVecNumeric
     * @brief Represents the element-wise multiplication of two `VecNumeric`
     * vectors.
     */
    class VecNumericTimesVecNumeric
        : public VecOperationSpecialization<VecNumericTimesVecNumeric>
    {
    public:
        VecNumericTimesVecNumeric(VecNumeric &a, VecNumeric &b) : a_(a), b_(b) {
            fatrop_dbg_assert(a.m() == b.m() && "Vector sizes must match");
        };
        Scalar operator()(const Index i) const { return a_(i) * b_(i); }
        Index m() const { return a_.m(); }
        friend VecNumericTimesVecNumeric operator*(VecNumeric &a, VecNumeric &b);
        VecNumeric& a() const { return a_; }
        VecNumeric& b() const { return b_; }
    private:
        VecNumeric &a_;
        VecNumeric &b_;
    };

    // operator overloading for VecNumeric specializations - blasfeo kernels
    template <>
    void VecNumeric::operator=(VecOperationSpecialization<VecNumericPlusVecNumeric>
                                   &&vecnumericplusvecnumeric)
    {
        auto vv = vecnumericplusvecnumeric.derived();
        fatrop_dbg_assert(this->m() == vv.a().m() && "Vector sizes must match");
        VecNumeric &a = vv.a();
        VecNumeric &b = vv.b();
        AXPY(m_, 1.0, &a.vec(), a.ai(), &b.vec(), b.ai(), &this->vec(), this->ai());
    }

    template <>
    void VecNumeric::operator=(
        VecOperationSpecialization<VecNumericTimesScalar> &&vecnumerictimesscalar)
    {
        auto vv = vecnumerictimesscalar.derived();
        fatrop_dbg_assert(this->m() == vv.a().m() && "Vector sizes must match");
        VecNumeric &a = vv.a();
        Scalar alpha = vv.alpha();
        VECCPSC(m_, alpha, &a.vec(), a.ai(), &this->vec(), this->ai());
    }

    template <>
    void VecNumeric::operator=(
        VecOperationSpecialization<VecNumericPlusVecNumericTimesScalar>
            &&vecnumericplusvecnumerictimesscalar)
    {
        auto vv = vecnumericplusvecnumerictimesscalar.derived();
        fatrop_dbg_assert(this->m() == vv.a().m() && "Vector sizes must match");
        VecNumeric &a = vv.a();
        VecNumeric &b = vv.b();
        Scalar alpha = vv.alpha();
        AXPY(m_, alpha, &b.vec(), b.ai(), &a.vec(), a.ai(), &this->vec(), this->ai());
    }

    template <>
    void VecNumeric::operator=(
        VecOperationSpecialization<VecNumericTimesScalarPlusVecNumericTimesScalar>
            &&vecnumerictimesscalarplusvecnumerictimesscalar)
    {
        auto vv = vecnumerictimesscalarplusvecnumerictimesscalar.derived();
        fatrop_dbg_assert(this->m() == vv.a().m() && "Vector sizes must match");
        VecNumeric &a = vv.a();
        Scalar alpha = vv.alpha();
        VecNumeric &b = vv.b();
        Scalar beta = vv.beta();
        AXPBY(m_, alpha, &a.vec(), a.ai(), beta, &b.vec(), b.ai(), &this->vec(),
              this->ai());
    }

    template <>
    void VecNumeric::operator=(VecOperationSpecialization<VecNumericTimesVecNumeric>
                                   &&vecnumerictimesscalar)
    {
        auto vv = vecnumerictimesscalar.derived();
        fatrop_dbg_assert(this->m() == vv.a().m() && "Vector sizes must match");
        VecNumeric &a = vv.a();
        VecNumeric &b = vv.b();
        VECMUL(m_, &b.vec(), b.ai(), &a.vec(), a.ai(), &this->vec(), this->ai());
    }

    /**
     * @brief Addition operator for VecNumeric and VecNumericTimesScalar.
     */
    VecNumericPlusVecNumericTimesScalar operator+(VecNumeric &a,
                                                  const VecNumericTimesScalar &b)
    {
        return VecNumericPlusVecNumericTimesScalar(a, b.a(), b.alpha());
    }

    /**
     * @brief Addition operator for VecNumericTimesScalar and VecNumeric.
     */
    VecNumericPlusVecNumericTimesScalar operator+(const VecNumericTimesScalar &a,
                                                  VecNumeric &b)
    {
        return b + a;
    }

    /**
     * @brief Addition operator for two VecNumericTimesScalar objects.
     */
    VecNumericTimesScalarPlusVecNumericTimesScalar operator+(
        const VecNumericTimesScalar &a, const VecNumericTimesScalar &b)
    {
        return VecNumericTimesScalarPlusVecNumericTimesScalar(a.a(), a.alpha(), b.a(),
                                                              b.alpha());
    }

    /**
     * @brief Element-wise multiplication operator for two VecNumeric objects.
     */
    VecNumericTimesVecNumeric operator*(VecNumeric &a, VecNumeric &b)
    {
        return VecNumericTimesVecNumeric(a, b);
    }

    /**
     * @brief Addition operator for two VecNumeric objects.
     */
    VecNumericPlusVecNumeric operator+(VecNumeric &a, VecNumeric &b)
    {
        return VecNumericPlusVecNumeric(a, b);
    }

    /**
     * @brief Scalar multiplication operator for VecNumeric.
     */
    VecNumericTimesScalar operator*(const Scalar alpha, VecNumeric &a)
    {
        return VecNumericTimesScalar(a, alpha);
    }

    /**
     * @brief Scalar multiplication operator for VecNumeric.
     */
    VecNumericTimesScalar operator*(VecNumeric &a, const Scalar alpha)
    {
        return alpha * a;
    }

    /**
     * @brief Unary minus operator for VecNumeric.
     */
    VecNumericTimesScalar operator-(VecNumeric &a) { return -1. * a; }

} // namespace fatrop

#endif // __fatrop_vector_specialization_hpp__
