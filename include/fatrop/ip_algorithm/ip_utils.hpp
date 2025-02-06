//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_utils__
#define __fatrop_ip_algorithm_ip_utils__
#include "fatrop/context/context.hpp"

namespace fatrop
{
    namespace internal
    {
        /**
         * @brief Compare two scalars with a relative tolerance.
         * 
         * This function compares two scalars 'a' and 'b' to determine if 'a' is less than or equal to 'b',
         * taking into account a relative tolerance based on the 'base' value.
         *
         * @param a The first scalar to compare.
         * @param b The second scalar to compare.
         * @param base The base value used to determine the relative tolerance.
         * @return bool True if 'a' is less than or equal to 'b' within the relative tolerance, false otherwise.
         */
        bool compare_le(const Scalar a, const Scalar b, const Scalar base);
    };
} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_utils__
