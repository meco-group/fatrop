#ifndef FATROPAUXINCLUDED
#define FATROPAUXINCLUDED
#include <vector>
/** \brief function to sum two integer vectors */
vector<int> operator+(const vector<int> a, const vector<int> b)
{
    vector<int> res = a;
    for (long unsigned int i = 0; i < a.size(); i++)
    {
        res.at(i) += b.at(i);
    }
    return res;
}
/** \brief function to sum integer vector and constant */
vector<int> operator+(const vector<int> a, const int b)
{
    vector<int> res = a;
    for (long unsigned int i = 0; i < a.size(); i++)
    {
        res.at(i) += b;
    }
    return res;
}
#endif //FATROPAUXINCLUDED