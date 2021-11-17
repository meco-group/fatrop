#ifndef FATROPAUXINCLUDED
#define FATROPAUXINCLUDED
#include <vector>
/** \brief function to sum two integer vectors */
vector<int> operator+(const vector<int> &a, const vector<int> &b)
{
    vector<int> res = a;
    for (long unsigned int i = 0; i < a.size(); i++)
    {
        res.at(i) += b.at(i);
    }
    return res;
}
/** \brief function to sum integer vector and constant */
vector<int> operator+(const vector<int> &a, const int b)
{
    vector<int> res = a;
    for (long unsigned int i = 0; i < a.size(); i++)
    {
        res.at(i) += b;
    }
    return res;
}
/** \brief function to cumulative sum integer vector, first element is zero */
vector<int> csum(const vector<int> &a)
{
    vector<int> res(a.size(), 0);
    for (long unsigned int i = 1; i < a.size(); i++)
    {
        res.at(i) = a.at(i) + res.at(i - 1);
    }
    return res;
}
int max(const vector<int> &a)
{
    int res = 0;
    for (long unsigned int i = 0; i < a.size(); i++)
    {
        int ai = a.at(i);
        res = abs(ai) > res ? ai : res;
    }
    return res;
}
#endif //FATROPAUXINCLUDED