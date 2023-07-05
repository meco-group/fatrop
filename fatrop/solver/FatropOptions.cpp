/*
 * Fatrop - A fast trajectory optimization solver
 * Copyright (C) 2022, 2023 Lander Vanroye, KU Leuven. All rights reserved.
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
#include "solver/FatropOptions.hpp"
using namespace fatrop;
using namespace std;

Option<bool>::Option(const string &name, const string &description, bool *value, bool default_value) : name_(name), description_(description), value(value), default_value_(default_value){};
void Option<bool>::set(const bool &new_value)
{
    *value = new_value;
}
void Option<bool>::set_default() const
{
    *value = default_value_;
}

template <typename T>
Option<T>::Option(const string &name, const string &description, T *value, T default_value, bool lower_bound_inclusive, T lower_bound, bool upper_bound_inclusive, T upper_bound) : name_(name), description_(description), value(value), default_value_(default_value), lower_bound_inclusive_(lower_bound_inclusive), lower_bound_(lower_bound), upper_bound_inclusive_(upper_bound_inclusive), upper_bound_(upper_bound){};
template <typename T>
Option<T> Option<T>::lower_bounded(const string &name, const string &description, T *value, T default_value, T lower_bound)
{
    return Option<T>(name, description, value, default_value, true, lower_bound, false, 0.0);
};
template <typename T>
Option<T> Option<T>::upper_bounded(const string &name, const string &description, T *value, T default_value, T upper_bound)
{
    return Option<T>(name, description, value, default_value, false, 0.0, true, upper_bound);
};
template <typename T>
Option<T> Option<T>::un_bounded(const string &name, const string &description, T *value, T default_value)
{
    return Option<T>(name, description, value, default_value, false, 0.0, false, 0.0);
};
template <typename T>
Option<T> Option<T>::box_bounded(const string &name, const string &description, T *value, T default_value, T lower_bound, T upper_bound)
{
    return Option<T>(name, description, value, default_value, true, lower_bound, true, upper_bound);
};
template <typename T>
void Option<T>::set_default() const
{
    *value = default_value_;
};
template <typename T>
void Option<T>::set(const T &new_value)
{
    // check if new value is in bounds
    if (lower_bound_inclusive_ && new_value < lower_bound_)
    {
        throw runtime_error("Option " + name_ + " is out of bounds");
    }
    if (upper_bound_inclusive_ && new_value > upper_bound_)
    {
        throw runtime_error("Option " + name_ + " is out of bounds");
    }
    *value = new_value;
};
FatropOptions::FatropOptions()
{
    register_option(IntegerOption::box_bounded("max_iter", "maximum number of iterations", &maxiter, 1000, 0, maxiter));
    register_option(NumericOption::lower_bounded("kappa_d", "kappa_d", &kappa_d, 1e-5, 0.0));
};
template <typename T>
void FatropOptions::set(const string &option_name, T value)
{
    if (numeric_options.find(option_name) != numeric_options.end())
    {
        if constexpr (std::is_floating_point<T>::value)
        {
            numeric_options[option_name].set(value);
        }
        else
        {
            throw std::runtime_error("Option " + option_name + " of type double");
        }
    }
    else if (integer_options.find(option_name) != integer_options.end())
    {
        if constexpr (std::is_integral<T>::value)
        {
            integer_options[option_name].set(value);
        }
        else if constexpr (std::is_floating_point<T>::value)
        {
            if ((int)value == value)
            {
                integer_options[option_name].set((int)value);
            }
            else
            {
                throw std::runtime_error("Option " + option_name + " of type int");
            }
        }
        else
        {
            throw std::runtime_error("Option " + option_name + " of type int");
        }
    }
    else if (boolean_options.find(option_name) != boolean_options.end())
    {
        if constexpr (std::is_same<T, bool>::value)
        {
            boolean_options[option_name].set(value);
        }
        else if constexpr (std::is_integral<T>::value || std::is_floating_point<T>::value)
        {
            if (value == 0)
            {
                boolean_options[option_name].set(false);
            }
            else if (value == 1)
            {
                boolean_options[option_name].set(true);
            }
            else
            {
                throw std::runtime_error("Option " + option_name + " of type bool can only convert integer value 0 or 1 to bool, got " + to_string(value));
            }
        }
        else
        {
            throw std::runtime_error("Option " + option_name + " of type bool");
        }
    }
    else
    {
        throw std::runtime_error("Option " + option_name + " not found");
    }
}

void FatropOptions::register_option(const NumericOption &option)
{
    numeric_options[option.name_] = option;
    option.set_default();
}
void FatropOptions::register_option(const IntegerOption &option)
{
    integer_options[option.name_] = option;
    option.set_default();
}
void FatropOptions::register_option(const BooleanOption &option)
{
    boolean_options[option.name_] = option;
    option.set_default();
}
auto operator<<(std::ostream &os, const FatropOptions &m) -> std::ostream &
{
    os << "Numeric options :" << std::endl;
    {
        for (auto const &x : m.numeric_options)
        {
            os << "   " << x.first << " : " << *x.second.value << std::endl;
        }
    }
    os << "Integer options :" << std::endl;
    {
        for (auto const &x : m.integer_options)
        {
            os << "   " << x.first << " : " << *x.second.value << std::endl;
        }
    }
    os << "Boolean options :" << std::endl;
    {
        for (auto const &x : m.boolean_options)
        {
            os << "   " << x.first << " : " << *x.second.value << std::endl;
        }
    }
    return os;
}
template class Option<fatrop_int>;
template class Option<double>;
template void FatropOptions::set<double>(const string &, double);
template void FatropOptions::set<fatrop_int>(const string &, int);
template void FatropOptions::set<bool>(const string &, bool);