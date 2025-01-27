//
// Copyright (C) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_filter_hpp__
#define __fatrop_ip_algorithm_ip_filter_hpp__
#include "fatrop/context/context.hpp"
#include <vector>

namespace fatrop
{
    struct IpFilterData
    {
        Scalar obj_value;
        Scalar constr_viol;
    };

    class IpFilter
    {
    public:
        IpFilter();
        void reset();
        void reserve(const Index max_size);
        bool is_acceptable(const IpFilterData &data) const;
        void augment(const IpFilterData &data);
        Index size() const;

    private:
        static bool a_dominates_b(const IpFilterData &a, const IpFilterData &b);
        std::vector<IpFilterData> data_;
    };
} // namespace fatrop

#endif //__fatrop_ip_algorithm_ip_algorithm_hpp__
