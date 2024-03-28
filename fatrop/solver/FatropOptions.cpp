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
#include "fatrop/solver/FatropOptions.hpp"
using namespace fatrop;
using namespace std;

template <typename T>
Option<T>::Option(const string &name, const string &description, T *value, T default_value) : name_(name), description_(description), value(value), default_value_(default_value){};
template <typename T>
void Option<T>::set(const T &new_value) const
{
    *value = new_value;
}
template <typename T>
void Option<T>::set_default() const
{
    *value = default_value_;
}

template <typename T>
NumberOption<T>::NumberOption(const string &name, const string &description, T *value, T default_value, bool lower_bound_inclusive, T lower_bound, bool upper_bound_inclusive, T upper_bound) : Option<T>(name, description, value, default_value), lower_bound_inclusive_(lower_bound_inclusive), lower_bound_(lower_bound), upper_bound_inclusive_(upper_bound_inclusive), upper_bound_(upper_bound){};
template <typename T>
NumberOption<T> NumberOption<T>::lower_bounded(const string &name, const string &description, T *value, T default_value, T lower_bound)
{
    return NumberOption<T>(name, description, value, default_value, true, lower_bound, false, 0.0);
};
template <typename T>
NumberOption<T> NumberOption<T>::upper_bounded(const string &name, const string &description, T *value, T default_value, T upper_bound)
{
    return NumberOption<T>(name, description, value, default_value, false, 0.0, true, upper_bound);
};
template <typename T>
NumberOption<T> NumberOption<T>::un_bounded(const string &name, const string &description, T *value, T default_value)
{
    return NumberOption<T>(name, description, value, default_value, false, 0.0, false, 0.0);
};
template <typename T>
NumberOption<T> NumberOption<T>::box_bounded(const string &name, const string &description, T *value, T default_value, T lower_bound, T upper_bound)
{
    return NumberOption<T>(name, description, value, default_value, true, lower_bound, true, upper_bound);
};
template <typename T>
void NumberOption<T>::set(const T &new_value) const
{
    // check if new value is in bounds
    if (lower_bound_inclusive_ && new_value < lower_bound_)
    {
        throw runtime_error("Option " + this->name_ + " is out of bounds");
    }
    if (upper_bound_inclusive_ && new_value > upper_bound_)
    {
        throw runtime_error("Option " + this->name_ + " is out of bounds");
    }
    *(this->value) = new_value;
};

FatropOptions::FatropOptions()
{
    register_option(IntegerOption::box_bounded("max_iter", "maximum number of iterations", &maxiter, 1000, 0, maxiter));
    register_option(DoubleOption::lower_bounded("kappa_d", "kappa_d", &kappa_d, 1e-5, 0.0));
};

bool FatropOptions::has_option(const std::string &option_name) const
{
    if (numeric_options.find(option_name) != numeric_options.end())
    {
        return true;
    }
    if(integer_options.find(option_name) != integer_options.end())
    {
        return true;
    }
    if (boolean_options.find(option_name) != boolean_options.end())
    {
        return true;
    }
    if (string_options.find(option_name) != string_options.end())
    {
        return true;
    }
    return false;
}
template <typename T>
void FatropOptions::set(const string &option_name, T value) const
{
    if (numeric_options.find(option_name) != numeric_options.end())
    {
        if constexpr (std::is_floating_point<T>::value)
        {
            numeric_options.at(option_name).set(value);
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
            integer_options.at(option_name).set(value);
        }
        else if constexpr (std::is_floating_point<T>::value)
        {
            if ((int)value == value)
            {
                integer_options.at(option_name).set((int)value);
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
            boolean_options.at(option_name).set(value);
        }
        else if constexpr (std::is_integral<T>::value || std::is_floating_point<T>::value)
        {
            if (value == 0)
            {
                boolean_options.at(option_name).set(false);
            }
            else if (value == 1)
            {
                boolean_options.at(option_name).set(true);
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
template <typename T>
void FatropOptions::prebuilt_set(const string &option_name, T value)
{
    if constexpr (std::is_integral<T>::value || std::is_floating_point<T>::value)
    {
        prebuilt_double[option_name] = value;
    }
    if constexpr (std::is_same<T, std::string>::value)
    {
        prebuilt_string[option_name] = value;
    }
}

void FatropOptions::register_option(const DoubleOption &option)
{
    numeric_options[option.name_] = option;
    option.set_default();
    // check if available in prebuilt options
    if (prebuilt_double.find(option.name_) != prebuilt_double.end())
    {
        option.set(prebuilt_double[option.name_]);
    }
}
void FatropOptions::register_option(const IntegerOption &option)
{
    integer_options[option.name_] = option;
    option.set_default();
    // check if available in prebuilt options
    if (prebuilt_double.find(option.name_) != prebuilt_double.end())
    {
        option.set(prebuilt_double[option.name_]);
    }
}
void FatropOptions::register_option(const BooleanOption &option)
{
    boolean_options[option.name_] = option;
    option.set_default();
    // check if available in prebuilt options
    if (prebuilt_double.find(option.name_) != prebuilt_double.end())
    {
        option.set(prebuilt_double[option.name_]);
    }
}

void FatropOptions::register_option(const StringOption &option)
{
    string_options[option.name_] = option;
    option.set_default();
    // check if available in prebuilt options
    if (prebuilt_string.find(option.name_) != prebuilt_string.end())
    {
        option.set(prebuilt_string[option.name_]);
    }
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
template class fatrop::NumberOption<fatrop_int>;
template class fatrop::NumberOption<double>;
template class fatrop::Option<std::string>;
template class fatrop::Option<bool>;
template void FatropOptions::set<double>(const string &, double) const;
template void FatropOptions::set<fatrop_int>(const string &, int) const;
template void FatropOptions::set<bool>(const string &, bool) const;
template void FatropOptions::set<std::string>(const string &, std::string) const;
template void FatropOptions::prebuilt_set<double>(const string &, double);
template void FatropOptions::prebuilt_set<int>(const string &, int);
template void FatropOptions::prebuilt_set<bool>(const string &, bool);
template void FatropOptions::prebuilt_set<std::string>(const string &, std::string);