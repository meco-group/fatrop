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
        const int size_a = a.size();
        FatropVector<T> res(size_a);
        res.at(0) = 0;
        for (int i = 1; i < size_a; i++)
        {
            res.at(i) = a.get(i - 1) + res.get(i - 1);
        }
        return res;
    }

    /** \brief returns index of max el of VecExpr */
    template <typename T, typename E>
    int max(const VecExpr<E, T> &a)
    {
        const int size_a = a.size();
        int res = 0;
        for (int i = 0; i < size_a; i++)
        {
            int ai = a.get(i);
            res = ai > res ? ai : res;
        }
        return res;
    }
};     // namespace fatrop
#endif // FATROPAUXINCLUDED