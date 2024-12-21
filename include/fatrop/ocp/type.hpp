//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ocp_type_hpp__
#define __fatrop_ocp_type_hpp__

/**
 * @file type.hpp
 * @brief Defines the OcpType structure for Optimal Control Problems.
 * 
 * This file contains the declaration of the OcpType structure,
 * which is used as a base for problem class-specific specializations
 */

namespace fatrop
{
    /**
     * @brief Structure representing the type of Optimal Control Problems.
     * 
     * This structure serves as a base for specializations of problem class-specific
     * elements such as Jacobians and Hessians. These specializations depend on the
     * problem's class structure and are crucial for efficient implementation of
     * optimization algorithms for different types of Optimal Control Problems.
     */
    struct OcpType
    {
    };
} // namespace fatrop

#endif //__fatrop_ocp_type_hpp__
