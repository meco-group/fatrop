//
// Copyright Lander Vanroye, KU Leuven
//

#ifndef __fatrop_linear_algebra_vector_specialization_hpp__

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
        VecNumericPlusVecNumeric(VecNumeric &a, VecNumeric &b) : a(a), b(b) {};
        Scalar operator[](const Index i) const { return a[i] + b[i]; }
        Index m() const { return a.m(); }
        VecNumeric &a;
        VecNumeric &b;
        friend VecNumericPlusVecNumeric operator+(VecNumeric &a, VecNumeric &b);
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
            : a(a), alpha(alpha) {};
        Scalar operator[](const Index i) const { return alpha * a[i]; }
        Index m() const { return a.m(); }
        VecNumeric &a;
        const Scalar alpha;
        friend VecNumericTimesScalar operator*(const Scalar alpha, VecNumeric &a);
        friend VecNumericTimesScalar operator*(VecNumeric &a, const Scalar alpha);
        friend VecNumericTimesScalar operator-(VecNumeric &a);
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
            : a(a), alpha(alpha), b(b), beta(beta) {};
        Scalar operator[](const Index i) const { return alpha * a[i] + beta * b[i]; }
        Index m() const { return a.m(); }
        VecNumeric &a;
        const Scalar alpha;
        VecNumeric &b;
        const Scalar beta;
        friend VecNumericTimesScalarPlusVecNumericTimesScalar operator+(
            const VecNumericTimesScalar &a, const VecNumericTimesScalar &b);
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
            : a(a), b(b), alpha(alpha) {};
        Scalar operator[](const Index i) const { return a[i] + alpha * b[i]; }
        Index m() const { return a.m(); }
        VecNumeric &a;
        VecNumeric &b;
        const Scalar alpha;
        friend VecNumericPlusVecNumericTimesScalar operator+(
            const VecNumericTimesScalar &a, VecNumeric &b);
        friend VecNumericPlusVecNumericTimesScalar operator+(
            VecNumeric &a, const VecNumericTimesScalar &b);
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
        VecNumericTimesVecNumeric(VecNumeric &a, VecNumeric &b) : a(a), b(b) {};
        Scalar operator[](const Index i) const { return a[i] * b[i]; }
        Index m() const { return a.m(); }
        VecNumeric &a;
        VecNumeric &b;
        friend VecNumericTimesVecNumeric operator*(VecNumeric &a, VecNumeric &b);
    };

    // operator overloading for VecNumeric specializations - blasfeo kernels
    template <>
    void VecNumeric::operator=(VecOperationSpecialization<VecNumericPlusVecNumeric>
                                   &&vecnumericplusvecnumeric)
    {
        auto vv = vecnumericplusvecnumeric.derived();
        VecNumeric &a = vv.a;
        VecNumeric &b = vv.b;
        AXPY(m_, 1.0, &a.vec(), a.ai(), &b.vec(), b.ai(), &this->vec(), this->ai());
    }

    template <>
    void VecNumeric::operator=(
        VecOperationSpecialization<VecNumericTimesScalar> &&vecnumerictimesscalar)
    {
        auto vv = vecnumerictimesscalar.derived();
        VecNumeric &a = vv.a;
        Scalar alpha = vv.alpha;
        VECCPSC(m_, alpha, &a.vec(), a.ai(), &this->vec(), this->ai());
    }

    template <>
    void VecNumeric::operator=(
        VecOperationSpecialization<VecNumericPlusVecNumericTimesScalar>
            &&vecnumericplusvecnumerictimesscalar)
    {
        auto vv = vecnumericplusvecnumerictimesscalar.derived();
        VecNumeric &a = vv.a;
        VecNumeric &b = vv.b;
        Scalar alpha = vv.alpha;
        AXPY(m_, alpha, &b.vec(), b.ai(), &a.vec(), a.ai(), &this->vec(), this->ai());
    }

    template <>
    void VecNumeric::operator=(
        VecOperationSpecialization<VecNumericTimesScalarPlusVecNumericTimesScalar>
            &&vecnumerictimesscalarplusvecnumerictimesscalar)
    {
        auto vv = vecnumerictimesscalarplusvecnumerictimesscalar.derived();
        VecNumeric &a = vv.a;
        Scalar alpha = vv.alpha;
        VecNumeric &b = vv.b;
        Scalar beta = vv.beta;
        AXPBY(m_, alpha, &a.vec(), a.ai(), beta, &b.vec(), b.ai(), &this->vec(),
              this->ai());
    }

    template <>
    void VecNumeric::operator=(VecOperationSpecialization<VecNumericTimesVecNumeric>
                                   &&vecnumerictimesscalar)
    {
        auto vv = vecnumerictimesscalar.derived();
        VecNumeric &a = vv.a;
        VecNumeric &b = vv.b;
        VECMUL(m_, &b.vec(), b.ai(), &a.vec(), a.ai(), &this->vec(), this->ai());
    }

    VecNumericPlusVecNumericTimesScalar operator+(VecNumeric &a,
                                                  const VecNumericTimesScalar &b)
    {
        return VecNumericPlusVecNumericTimesScalar(a, b.a, b.alpha);
    }

    VecNumericPlusVecNumericTimesScalar operator+(const VecNumericTimesScalar &a,
                                                  VecNumeric &b)
    {
        return b + a;
    }

    VecNumericTimesScalarPlusVecNumericTimesScalar operator+(
        const VecNumericTimesScalar &a, const VecNumericTimesScalar &b)
    {
        return VecNumericTimesScalarPlusVecNumericTimesScalar(a.a, a.alpha, b.a,
                                                              b.alpha);
    }

    VecNumericTimesVecNumeric operator*(VecNumeric &a, VecNumeric &b)
    {
        return VecNumericTimesVecNumeric(a, b);
    }

    VecNumericPlusVecNumeric operator+(VecNumeric &a, VecNumeric &b)
    {
        return VecNumericPlusVecNumeric(a, b);
    }

    VecNumericTimesScalar operator*(const Scalar alpha, VecNumeric &a)
    {
        return VecNumericTimesScalar(a, alpha);
    }

    VecNumericTimesScalar operator*(VecNumeric &a, const Scalar alpha)
    {
        return alpha * a;
    }

    VecNumericTimesScalar operator-(VecNumeric &a) { return -1. * a; }

} // namespace fatrop

#endif // __fatrop_vector_specialization_hpp__