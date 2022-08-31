#ifndef FILTERINCLUDED
#define FILTERINCLUDED
#include "vector"
#include "aux/SmartPtr.hpp"
using namespace std;
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
        void Augment(const FilterData &filterdata);
        /** \brief check if fdin0 is dominated by fdin1 */
        inline bool IsDominated(const FilterData &fdin0, const FilterData &fdin1) const;
        bool IsAcceptable(const FilterData &fdin) const;
        void Reset();

    private:
        vector<FilterData> filterdata_;
    };
} // namespace fatrop

#endif // FILTERINCLUDED