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
        virtual LinsolReturnFlag compute_search_dir() = 0;
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
        VecRealAllocated Di_;
        VecRealAllocated Ds_;
        VecRealAllocated Deq_;
    };
} // namespace fatrop

#endif