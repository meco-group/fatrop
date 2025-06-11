//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_common_fwd_hpp__
#define __fatrop_common_fwd_hpp__


namespace fatrop
{
    // Forward declarations for exception.hpp
    // class FatropException;

    // Forward declarations for options.hpp
    class OptionSetterBase;
    template <typename OptionType, typename AlgoType> class OptionSetter;
    class OptionRegistry;

    // Forward declarations for printing.hpp
    class OutputStreamManager;

    class Timer;

    // Forward declarations for timing.hpp
    // No forward declarations needed as it only includes external headers and defines macros

} // namespace fatrop

#endif // __fatrop_common_fwd_hpp__
