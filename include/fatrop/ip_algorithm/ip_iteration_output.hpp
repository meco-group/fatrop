//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_iteration_output_hpp__
#define __fatrop_ip_iteration_output_hpp__
#include "fatrop/context/context.hpp"
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include <memory>

namespace fatrop
{
    class IpIterationOutputBase
    {
    public:
        virtual ~IpIterationOutputBase() = default;
        virtual void print_header() = 0;
        virtual void output_current_iteration() = 0;
    };

    template<typename ProblemType>
    class IpIterationOutput : public IpIterationOutputBase
    {
    public:
        typedef std::shared_ptr<IpData<ProblemType>> IpDataSp;
        typedef IpIterate<ProblemType> IpIterateType;

        IpIterationOutput(const IpDataSp &ipdata);
        void print_header() override;
        void output_current_iteration() override;

    private:
        void print_iteration(Index iter, Scalar objective, Scalar inf_pr, Scalar inf_du,
                             Scalar lg_mu, Scalar d_norm, Scalar rg, Scalar alpha_du,
                             Scalar alpha_pr, Index ls);
        IpDataSp ipdata_;
    };

} // namespace fatrop

#endif // __fatrop_ip_iteration_output_hpp__
