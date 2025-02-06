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
    enum class IpConvergenceStatus
    {
        Continue,
        Converged,
        ConvergedToAcceptablePoint,
        MaxIterExceeded, 
        Unassigned
    };

    class IpConvergenceCheckBase
    {
    public:
        virtual IpConvergenceStatus check_converged() = 0;
        virtual bool check_acceptable() const = 0;
        virtual void reset() = 0;

    protected:
        virtual ~IpConvergenceCheckBase() = default;
    };

    template <typename ProblemType> class IpConvergenceCheck : public IpConvergenceCheckBase
    {
        typedef std::shared_ptr<IpData<ProblemType>> IpDataSp;

    public:
        IpConvergenceCheck(const IpDataSp &ipdata);
        IpConvergenceStatus check_converged() override;
        void reset() override;
        bool check_acceptable() const override;

    private:
        IpDataSp ipdata_;
        Scalar tol_acceptable_ = 1e-6;
        Index acceptable_counter_ = 0;
        Index acceptable_iter_ = 15;
        Index max_iter_ = 1000;
    };

} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_convergence_check_hpp__
