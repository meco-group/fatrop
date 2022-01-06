#ifndef FATROPAUXINCLUDED
#define FATROPAUXINCLUDED
#include <vector>
#include "FatropVector.hpp"
namespace fatrop
{
    /** \brief function to cumulative sum integer vector, first element is zero */
    template <typename T, typename E>
    FatropVector<T> offsets(const VecExpr<E, T> &a)
    {
        FatropVector<T> res(a.size());
        res.at(0) = 0;
        for (int i = 1; i < a.size(); i++)
        {
            res.at(i) = a.getEl(i - 1) + res.getEl(i - 1);
        }
        return res;
    }

    /** \brief returns index of max el of VecExpr*/
    template <typename T, typename E>
    int max(const VecExpr<E, T> &a)
    {
        int res = 0;
        for (int i = 0; i < a.size(); i++)
        {
            int ai = a.getEl(i);
            res = ai > res ? ai : res;
        }
        return res;
    }
};     // namespace fatrop
#endif // FATROPAUXINCLUDED