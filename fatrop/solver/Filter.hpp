#ifndef FILTERINCLUDED
#define FILTERINCLUDED
#include "vector"
namespace fatrop
{
    struct FilterData
    {
        FilterData(){};
        FilterData(const int iteration, const double phi, const double theta):iteration(iteration),phi(phi),theta(theta){};
        const int iteration = 0;
        /** \brief barrier function value */
        const double phi = 0.0;
        /** \brief constraint violation value */
        const double theta = 0.0;
    };
    class Filter
    {
    public:
        Filter(const int size);
        void augment(const FilterData &filterdata);
        /** \brief check if fdin0 is dominated by fdin1 */
        inline bool a_dominmates_b(const FilterData &fdin0, const FilterData &fdin1) const;
        bool is_acceptable(const FilterData &fdin) const;
        void reset();
        int size() const
        {
            return filterdata_.size();
        }

    private:
        std::vector<FilterData> filterdata_;
    };
} // namespace fatrop

#endif // FILTERINCLUDED