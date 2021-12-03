#ifndef FATROPAUXINCLUDED
#define FATROPAUXINCLUDED
#include <vector>
namespace fatrop
{
    // /** \brief function to sum two integer vectors */
    // vector<int> operator+(const vector<int> &a, const vector<int> &b)
    // {
    //     vector<int> res = a;
    //     for (long unsigned int i = 0; i < a.size(); i++)
    //     {
    //         res.at(i) += b.at(i);
    //     }
    //     return res;
    // }
    // /** \brief function to sum integer vector and constant */
    // vector<int> operator+(const vector<int> &a, const int b)
    // {
    //     vector<int> res = a;
    //     for (long unsigned int i = 0; i < a.size(); i++)
    //     {
    //         res.at(i) += b;
    //     }
    //     return res;
    // }
    /** \brief function to cumulative sum integer vector, first element is zero */

    template <typename T, typename E>
    FatropVector<T> offsets(const VecExpr<E,T> &a)
    {
        FatropVector<T> res(a.size());
        for (int i = 1; i < a.size(); i++)
        {
            res.at(i) = a.getEl(i) + res.getEl(i - 1);
        }
        return res;
    }

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
#endif //FATROPAUXINCLUDED