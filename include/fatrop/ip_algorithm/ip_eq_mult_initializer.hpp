//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_eq_mult_initializer_hpp__
#define __fatrop_ip_algorithm_ip_eq_mult_initializer_hpp__
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"

#include "fatrop/nlp/fwd.hpp"
#include <memory>

namespace fatrop
{
    class IpEqMultInitializerBase
    {
    public:
        virtual void initialize_eq_mult() = 0;
        virtual void reset() = 0;
    protected:
        virtual ~IpEqMultInitializerBase() = default;
    };

    template <typename ProblemType> class IpEqMultInitializer : public IpEqMultInitializerBase
    {
        typedef std::shared_ptr<PdSolverOrig<ProblemType>> PdSolverSp;
        typedef std::shared_ptr<IpData<ProblemType>> IpDataSp;
        typedef IpIterate<ProblemType> IpIterateType;
        typedef std::shared_ptr<IpIterateType> IpIterateSp;

    public:
        IpEqMultInitializer(const IpDataSp &ipdata, const PdSolverSp &linear_solver);
        void initialize_eq_mult() override;
        void reset() override;

    private:
        IpDataSp ipdata_;
        PdSolverSp linear_solver_;
        VecRealAllocated rhs_x_;
        VecRealAllocated rhs_s_;
        VecRealAllocated rhs_g_;
        VecRealAllocated rhs_cl_;
        VecRealAllocated rhs_cu_;
        VecRealAllocated Dx_;
        VecRealAllocated Ds_;
        VecRealAllocated Deq_;
        VecRealAllocated dummy_s_;
        VecRealAllocated dummy_z_;
        Scalar lam_max_;
    };
} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_eq_mult_initializer_hpp__
