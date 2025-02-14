//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_resto_phase_hpp__
#define __fatrop_ip_resto_phase_hpp__

#include "fatrop/common/fwd.hpp"
namespace fatrop
{
    class IpRestoPhaseBase
    {
        public:
        virtual bool perform_restoration() = 0;
        virtual void register_options(OptionRegistry &registry) = 0;
        virtual void reset() = 0;
        virtual ~IpRestoPhaseBase() = default;
    };
} // namespace fatrop

#endif // __fatrop_ip_resto_phase_hpp__
