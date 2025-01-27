#include "fatrop/ip_algorithm/ip_filter.hpp"

namespace fatrop
{
    IpFilter::IpFilter() {}

    void IpFilter::reset() { data_.resize(0); }

    void IpFilter::reserve(const Index max_size) { data_.reserve(max_size + 1); }

    bool IpFilter::is_acceptable(const IpFilterData &data) const
    {
        for (const IpFilterData &d : data_)
        {
            if (a_dominates_b(data, d))
                return false;
        }
        return true;
    }
    void IpFilter::augment(const IpFilterData &data) { data_.push_back(data); }

    Index IpFilter::size() const { return data_.size(); }

    bool IpFilter::a_dominates_b(const IpFilterData &a, const IpFilterData &b)
    {
        return (a.obj_value > b.obj_value && a.constr_viol > b.constr_viol);
    }
} // namespace fatrop
