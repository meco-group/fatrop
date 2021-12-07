/**
 * @file FatropVector.hpp
 * @author Lander Vanroye
 * @brief This file containts the FatropVector, which is derived from the std::vector but also allows templated expressions (implemented with CRTP - static polymorphism) such as vector-vector sum and vector-scalar sum. https://en.wikipedia.org/wiki/Expression_templates
 * @version 0.1
 * @date 2021-12-03
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef FATROPVECTORINCLUDED
#define FATROPVECTORINCLUDED
#include <vector>
#include <assert.h>
#include <utility>
using namespace std;
namespace fatrop
{
    // static polymorphism using CRTP
    template <typename E, typename T>
    class VecExpr
    {
    public:
        T getEl(const int ai) const { return static_cast<E const &>(*this).getEl(ai); };
        int size() const { return static_cast<E const &>(*this).size(); };
    };

    template <typename T, typename E1, typename E2>
    class VecSum : public VecExpr<VecSum<T, E1, E2>, T>
    {
    public:
        VecSum(const VecExpr<E1, T> &expr1, const VecExpr<E2, T> &expr2) : expr1_(expr1), expr2_(expr2)
        {
            assert(expr1.size() == expr2.size());
        };
        T getEl(const int ai) const { return expr1_.getEl(ai) + expr2_.getEl(ai); };
        int size() const { return expr1_.size(); };

    private:
        const VecExpr<E1, T> &expr1_;
        const VecExpr<E2, T> &expr2_;
    };
    template <typename T, typename E1>
    class VecScalarSum : public VecExpr<VecScalarSum<T, E1>, T>
    {
    public:
        VecScalarSum(const VecExpr<E1, T> &expr, const int scalar) : expr_(expr), scalar_(scalar){};
        T getEl(const int ai) const { return expr_.getEl(ai) + scalar_; };
        int size() const { return expr_.size(); };

    public:
        const VecExpr<E1, T> &expr_;
        T scalar_;
    };

    template <typename T>
    class FatropVector : public vector<T>, public VecExpr<FatropVector<T>, T>
    {
    public:
        FatropVector() : vector<T>(){};
        FatropVector(const int size) : vector<T>(size){};
        template <typename E>
        FatropVector(const VecExpr<E, T> &vecexpr) : vector<T>(vecexpr.size())
        {
            // todo: vector is first initialized, initialize with iterator
            for (int i = 0; i < vecexpr.size(); i++)
            {
                vector<T>::at(i) = vecexpr.getEl(i);
            }
        }
        FatropVector(const vector<T> &vec) : vector<T>(vec){};
        FatropVector(vector<T> &&vec) : vector<T>(move(vec)){};
        T getEl(const int ai) const { return vector<T>::at(ai); };
        int size() const { return vector<T>::size(); };
    };
    template <typename T, typename E1, typename E2>
    VecSum<T, E1, E2> operator+(const VecExpr<E1, T> &expr1, const VecExpr<E2, T> &expr2)
    {
        return VecSum<T, E1, E2>(expr1, expr2);
    }

    template <typename T, typename E1>
    VecScalarSum<T, E1> operator+(const VecExpr<E1, T> &expr1, const int scalar)
    {
        return VecScalarSum<T, E1>(expr1, scalar);
    }

    template <typename T, typename E1>
    VecScalarSum<T, E1> operator+(const int scalar, const VecExpr<E1, T> &expr1)
    {
        return VecScalarSum<T, E1>(expr1, scalar);
    }
} // namespace fatrop

#endif // FATROPVECTORINCLUDED
