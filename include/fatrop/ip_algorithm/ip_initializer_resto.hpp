//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_initializer_resto_hpp__
#define __fatrop_ip_algorithm_ip_initializer_resto_hpp__
#include "fatrop/common/fwd.hpp"
#include "fatrop/context/context.hpp"
#include "fatrop/ip_algorithm/fwd.hpp"
#include "fatrop/ip_algorithm/ip_initializer.hpp"
#include <memory>
namespace fatrop
{
    template <typename ProblemType> class IpInitializerResto : public IpInitializerBase
    {
        typedef IpData<ProblemType> IpDataType;
        typedef std::shared_ptr<IpDataType> IpDataSp;
        typedef IpIterate<ProblemType> IpIterateType;
        typedef ProblemInfo<ProblemType> InfoType;

    public:
        IpInitializerResto(const IpDataSp &data_orig, const IpDataSp &data_resto)
            : data_orig_(data_orig), data_resto_(data_resto)
        {
        }
        virtual void initialize() override;
        virtual void reset() override;
        void set_rho(const Scalar& rho) { rho_ = rho; }
        // Register options
        void register_options(OptionRegistry &registry) override;

    private:
        IpDataSp data_orig_;
        IpDataSp data_resto_;
        Scalar rho_ = 1000.;
    };

} // namespace fatrop

#endif // __fatrop_ip_algorithm_ip_initializer__resto_hpp__
