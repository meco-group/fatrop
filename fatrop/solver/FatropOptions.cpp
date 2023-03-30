#include "solver/FatropOptions.hpp"
using namespace fatrop;

Option<bool>::Option(const string &name, const string &description, bool *value, bool default_value) : name_(name), description_(description), value(value), default_value_(default_value){};
void Option<bool>::set(const bool &new_value)
{
    *value = new_value;
}
void Option<bool>::SetDefault() const
{
    *value = default_value_;
}

template <typename T>
Option<T>::Option(const string &name, const string &description, T *value, T default_value, bool lower_bound_inclusive, T lower_bound, bool upper_bound_inclusive, T upper_bound) : name_(name), description_(description), value(value), default_value_(default_value), lower_bound_inclusive_(lower_bound_inclusive), lower_bound_(lower_bound), upper_bound_inclusive_(upper_bound_inclusive), upper_bound_(upper_bound){};
template <typename T>
Option<T> Option<T>::LowerBounded(const string &name, const string &description, T *value, T default_value, T lower_bound)
{
    return Option<T>(name, description, value, default_value, true, lower_bound, false, 0.0);
};
template <typename T>
Option<T> Option<T>::UpperBounded(const string &name, const string &description, T *value, T default_value, T upper_bound)
{
    return Option<T>(name, description, value, default_value, false, 0.0, true, upper_bound);
};
template <typename T>
Option<T> Option<T>::UnBounded(const string &name, const string &description, T *value, T default_value)
{
    return Option<T>(name, description, value, default_value, false, 0.0, false, 0.0);
};
template <typename T>
Option<T> Option<T>::BoxBounded(const string &name, const string &description, T *value, T default_value, T lower_bound, T upper_bound)
{
    return Option<T>(name, description, value, default_value, true, lower_bound, true, upper_bound);
};
template <typename T>
void Option<T>::SetDefault() const
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
    RegisterOption(IntegerOption::BoxBounded("max_iter", "maximum number of iterations", &maxiter, 1000, 0, maxiter));
    RegisterOption(NumericOption::LowerBounded("kappa_d", "kappa_d", &kappa_d, 1e-5, 0.0));
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

void FatropOptions::RegisterOption(const NumericOption &option)
{
    numeric_options[option.name_] = option;
    option.SetDefault();
}
void FatropOptions::RegisterOption(const IntegerOption &option)
{
    integer_options[option.name_] = option;
    option.SetDefault();
}
void FatropOptions::RegisterOption(const BooleanOption &option)
{
    boolean_options[option.name_] = option;
    option.SetDefault();
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
template class Option<int>;
template class Option<double>;
template void FatropOptions::set<double>(const string &, double);
template void FatropOptions::set<int>(const string &, int);
template void FatropOptions::set<bool>(const string &, bool);