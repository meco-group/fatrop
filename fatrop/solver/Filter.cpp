#include "solver/Filter.hpp"
using namespace fatrop;
using namespace std;
Filter::Filter(const int size)
{
    filterdata_.reserve(size + 1);
}
void Filter::Augment(const FilterData &filterdata)
{
    filterdata_.push_back(filterdata);
}
inline bool Filter::IsDominated(const FilterData &fdin0, const FilterData &fdin1) const
{
    // worse barrier filter and constraint violation -> dominated
    if (fdin0.phi > fdin1.phi && fdin0.theta > fdin1.theta)
    {
        return true;
    }
    return false;
}
bool Filter::IsAcceptable(const FilterData &fdin) const
{
    // run over filterdata_ elements
    for (vector<double>::size_type k = 0; k < filterdata_.size(); k++)
    {
        if (IsDominated(fdin, filterdata_.at(k)))
        {
            return false;
        };
    }
    // fdin is not dominated by one of the filterdata_ elements -> acceptable to filter
    return true;
}
void Filter::Reset()
{
    filterdata_.resize(0);
}