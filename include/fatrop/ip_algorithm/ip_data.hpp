//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_data_hpp__
#define __fatrop_ip_algorithm_ip_data_hpp__
#include "fatrop/ip_algorithm/ip_iterate.hpp" 
#include "fatrop/context/context.hpp"
namespace fatrop
{
    template <typename ProblemType>
    struct IpData
    {
        IpIterate<ProblemType> iterates[3];
        IpIterate& it_curr = iterates[0];
        IpIterate& it_trial = iterates[1];
        IpIterate& it_backup = iterates[2];
    };

} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_data_hpp__