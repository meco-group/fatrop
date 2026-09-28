// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_IP_WORKSPACE_HPP
#define FATROP_IP_WORKSPACE_HPP

#include "fatrop/context/context.hpp"

namespace fatrop
{
    /// Storage budgets supplied before the algorithm's components are constructed.
    /// Runtime history options may use any size up to the reserved capacity.
    struct IpWorkspace
    {
        Index hessian_history_capacity = 6;
        Index rhs_capacity = 12;
        Index max_iterations_capacity = 1000;
    };
}
#endif
