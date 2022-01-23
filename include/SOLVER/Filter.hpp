#ifndef FILTERINCLUDED
#define FILTERINCLUDED
#include "vector"
using namespace std;
namespace fatrop
{
    struct FilterData
    {
        /** \brief barrier function value */
        double phi = 0.0;
        /** \brief constraint violation value */
        double theta = 0.0;
    };
    class Filter
    {
    public:
        Filter(const int size)
        {
            filterdata_.reserve(size + 1);
        }
        void Augment(const FilterData &filterdata)
        {
            filterdata_.push_back(filterdata);
        }
        /** \brief check if fdin0 is dominated by fdin1 */
        inline bool IsDominated(const FilterData &fdin0, const FilterData &fdin1) const
        {
            // worse barrier filter and constraint violation -> dominated
            if (fdin0.phi > fdin1.phi && fdin0.theta > fdin1.theta)
            {
                return true;
            }
            return false;
        }
        bool IsAcceptable(const FilterData &fdin)
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
        void Reset()
        {
            filterdata_.resize(0);
        }

    private:
        vector<FilterData> filterdata_;
    };
} // namespace fatrop

#endif // FILTERINCLUDED