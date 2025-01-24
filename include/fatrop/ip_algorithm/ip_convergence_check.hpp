//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_convergence_check_hpp__
#define __fatrop_ip_algorithm_ip_convergence_check_hpp__
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/linear_algebra/fwd.hpp"
#include "fatrop/linear_algebra/vector.hpp"
#include <memory>
namespace fatrop
{
    class IpConvergenceCheckBase
    {
    public:
        virtual bool converged() const = 0;

    protected:
        virtual ~IpConvergenceCheckBase() = default;
    };

    template <typename ProblemType> class IpConvergenceCheck : public IpConvergenceCheckBase
    {
        typedef std::shared_ptr<IpData<ProblemType>> IpDataSp;
    public:
        IpConvergenceCheck(const IpDataSp &ipdata);
        bool converged() const override;

    private:
        IpDataSp ipdata_;
        Scalar tol_ = 1e-8;
    };

} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_convergence_check_hpp__