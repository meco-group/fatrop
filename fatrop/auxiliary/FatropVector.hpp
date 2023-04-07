/*
 * Fatrop - A fast trajectory optimization solver
 * Copyright (C) 2022, 2023 Lander Vanroye <lander.vanroye@kuleuven.be>
 *
 * This file is part of Fatrop.
 *
 * Fatrop is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fatrop is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Fatrop.  If not, see <http://www.gnu.org/licenses/>. */
#ifndef FATROP_VECTOR_INCLUDED
#define FATROP_VECTOR_INCLUDED
#include <vector>
#include <assert.h>
#include <utility>
#include <functional>
namespace fatrop
{
    // static polymorphism using CRTP
    template <typename E, typename T>
    class VecExpr
    {
    public:
        T get(const int ai) const { return static_cast<E const &>(*this).get(ai); };
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
        T get(const int ai) const { return expr1_.get(ai) + expr2_.get(ai); };
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
        T get(const int ai) const { return expr_.get(ai) + scalar_; };
        int size() const { return expr_.size(); };

    public:
        const VecExpr<E1, T> &expr_;
        T scalar_;
    };
    template <typename T, typename E1>
    class VecRotate : public VecExpr<VecRotate<T, E1>, T>
    {
    public:
        VecRotate(const VecExpr<E1, T> &expr, const int shift) : expr_(expr), shift_(shift){};
        T get(const int ai) const { return expr_.get((ai + shift_ + (((ai + shift_) / size() + 1) * size())) % size()); };
        int size() const { return expr_.size(); };

    public:
        const VecExpr<E1, T> &expr_;
        int shift_;
    };

    template <typename T>
    class FatropVector : public std::vector<T>, public VecExpr<FatropVector<T>, T>
    {
    public:
        FatropVector() : std::vector<T>(){};
        FatropVector(const int size) : std::vector<T>(size){};
        template <typename E>
        FatropVector(const VecExpr<E, T> &vecexpr) : std::vector<T>(vecexpr.size())
        {
            // todo: vector is first initialized, initialize with iterator
            for (int i = 0; i < vecexpr.size(); i++)
            {
               std:: vector<T>::at(i) = vecexpr.get(i);
            }
        }
        FatropVector(const std::vector<T> &vec) : std::vector<T>(vec){};
        FatropVector(std::vector<T> &&vec) : std::vector<T>(move(vec)){};
        T get(const int ai) const { return std::vector<T>::at(ai); };
        int size() const { return std::vector<T>::size(); };
        operator T *() { return this->data(); };
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
    template <typename T, typename E1>
    T sum(const VecExpr<E1, T> &expr)
    {
        T res = 0;
        for (int i = 0; i < expr.size(); i++)
        {
            res += expr.get(i);
        }
        return res;
    }
    template <typename T, typename E1>
    VecRotate<T, E1> rotate(const VecExpr<E1, T> &expr, const int shift)
    {
        return VecRotate<T, E1>(expr, shift);
    };
    template <typename T>
    std::vector<T> TransformRange(const int begin, const int end, const std::function<T(int)>& func)
    {
        int size = end - begin;
        std::vector<T> res(size);
        for (int i = 0; i < size; i++)
        {
            res.at(i) = func(begin + i);
        }
        return res;
    }
} // namespace fatrop

#endif // FATROP_VECTOR_INCLUDED
