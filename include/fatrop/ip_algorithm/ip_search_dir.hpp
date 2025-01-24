//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_search_dir_hpp__
#define __fatrop_ip_algorithm_ip_search_dir_hpp__
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/linear_solver_return_flags.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include <memory>

namespace fatrop
{
    class IpSearchDirBase
    {
    public:
        virtual void reset() = 0;
        virtual LinsolReturnFlag compute_search_dir() = 0;
    protected:
        virtual ~IpSearchDirBase() = default;
    };

    template <typename ProblemType, typename LinearSystemType, typename LinearSolverDerived>
    class IpSearchDirImpl : public IpSearchDirBase
    {
        typedef IpIterate<ProblemType> IpIterateType;
        typedef std::shared_ptr<IpIterateType> IpIterateSp;
        typedef IpData<ProblemType> IpDataType;
        typedef std::shared_ptr<IpDataType> IpDataSp;
        typedef LinearSolver<LinearSolverDerived, LinearSystemType> LinearSolverType;
        typedef std::shared_ptr<LinearSolverType> LinearSolverSp;

    public:
        IpSearchDirImpl(const IpDataSp &ipdata, const LinearSolverSp &linear_solver);
        void reset() override;
        LinsolReturnFlag compute_search_dir() override;

    private:
        IpDataSp ipdata_;
        LinearSolverSp linear_solver_;
        VecRealAllocated rhs_x_;
        VecRealAllocated rhs_s_;
        VecRealAllocated rhs_g_;
        VecRealAllocated rhs_cl_;
        VecRealAllocated rhs_cu_;
        VecRealAllocated Dx_;
        VecRealAllocated Ds_;
        VecRealAllocated Deq_;
        Scalar delta_w_last_ = 0.;
        Scalar delta_w0_ = 1e-4;
        Scalar delta_wmin_ = 1e-20;
        Scalar kappa_wmin_ = 1. / 3.;
        Scalar kappa_wplus_ = 8.;
        Scalar kappa_wplusem_ = 100.;
        Scalar kappa_c_ = 0.25;
        Scalar delta_c_stripe_ = 1e-6;
    };
} // namespace fatrop

#endif