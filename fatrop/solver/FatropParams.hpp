#ifndef FATROPPARAMSINCLUDED
#define FATROPPARAMSINCLUDED
#include "aux/SmartPtr.hpp"
#include <string>
#include <map>
#include <string>
#include <type_traits>
using namespace std;
namespace fatrop
{
    template <typename T>
    struct Option
    {
    public:
        // NumericOption operator=(const NumericOption &other) = default;
        Option(){};
        Option(const string &name, const string &description, T *value, T default_value, bool lower_bound_inclusive, T lower_bound, bool upper_bound_inclusive, T upper_bound) : name_(name), description_(description), value(value), default_value_(default_value), lower_bound_inclusive_(lower_bound_inclusive), lower_bound_(lower_bound), upper_bound_inclusive_(upper_bound_inclusive), upper_bound_(upper_bound){};
        static Option<T> LowerBounded(const string &name, const string &description, T *value, T default_value, T lower_bound)
        {
            return Option<T>(name, description, value, default_value, true, lower_bound, false, 0.0);
        };
        static Option<T> UpperBounded(const string &name, const string &description, T *value, T default_value, T upper_bound)
        {
            return Option<T>(name, description, value, default_value, false, 0.0, true, upper_bound);
        };
        static Option<T> UnBounded(const string &name, const string &description, T *value, T default_value)
        {
            return Option<T>(name, description, value, default_value, false, 0.0, false, 0.0);
        };
        static Option<T> BoxBounded(const string &name, const string &description, T *value, T default_value, T lower_bound, T upper_bound)
        {
            return Option<T>(name, description, value, default_value, true, lower_bound, true, upper_bound);
        };
        void SetDefault() const
        {
            *value = default_value_;
        };
        void Set(const T &new_value)
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
        string name_;
        string description_;
        T *value = NULL;
        T default_value_;
        bool lower_bound_inclusive_;
        T lower_bound_;
        bool upper_bound_inclusive_;
        T upper_bound_;
    };
    // template <>
    // struct Option<string>
    // {
    //     Option<string>(const string& name, const string& descrription, )
    //     string name_;
    //     // void Set(const string& new_value)
    //     // {
    //     //     *value = new_value;
    //     // }
    // };
    template <>
    struct Option<bool>
    {
        Option<bool>(){};
        Option<bool>(const string &name, const string &description, bool *value, bool default_value) : name_(name), description_(description), value(value), default_value_(default_value){};
        void Set(const bool& new_value)
        {
            *value = new_value;
        }
        void SetDefault() const
        {
            *value = default_value_;
        }
        string name_;
        string description_;
        bool *value = NULL;
        bool default_value_;
    };
    // define Numeric option as Option<double>
    typedef Option<double> NumericOption;
    typedef Option<int> IntegerOption;
    typedef Option<string> StringOption;
    typedef Option<bool> BooleanOption;

    class FatropOptions
    {
    public:
        FatropOptions()
        {
            RegisterOption(IntegerOption::BoxBounded("max_iter", "maximum number of iterations", &maxiter, 1000, 0, maxiter));
            RegisterOption(NumericOption::LowerBounded("kappa_d", "kappa_d", &kappa_d, 1e-5, 0.0));
        };
        int maxiter = 1000; // TODO this value cannot be changed to a value larger than the one used for building the solver
        double kappa_d = 1e-5;
        template <typename T>
        void SetOption(const string &option_name, T value)
        {
            if (numeric_options.find(option_name) != numeric_options.end())
            {
                if constexpr (std::is_floating_point<T>::value)
                {
                    numeric_options[option_name].Set(value);
                }
                else
                {
                    throw std::runtime_error("Option " + option_name + " not of type double");
                }
            }
            else if (integer_options.find(option_name) != integer_options.end())
            {
                if constexpr (std::is_integral<T>::value)
                {
                    integer_options[option_name].Set(value);
                }
                else
                {
                    throw std::runtime_error("Option " + option_name + " not of type int");
                }
            }
            // else if (string_options.find(option_name) != string_options.end())
            // {
            //     if const_expr(std::is_same<T, string>::value())
            //     {
            //         string_options[option_name].Set(value);
            //     }
            //     else
            //     {
            //         throw std::runtime_error("Option " + option_name + " not of type string");
            //     }
            // }
            else if (boolean_options.find(option_name) != boolean_options.end())
            {
                if constexpr (std::is_same<T, bool>::value)
                {
                    boolean_options[option_name].Set(value);
                }
            }
            else
            {
                throw std::runtime_error("Option " + option_name + " not found");
            }
        }

    public:
        void RegisterOption(const NumericOption &option)
        {
            numeric_options[option.name_] = option;
            option.SetDefault();
        }
        void RegisterOption(const IntegerOption &option)
        {
            integer_options[option.name_] = option;
            option.SetDefault();
        }
        // void RegisterStringOption(const StringOption &option)
        // {
        //     string_options[option.name_] = option;
        // }
        void RegisterOption(const BooleanOption &option)
        {
            boolean_options[option.name_] = option;
            option.SetDefault();
        }
        map<string, NumericOption> numeric_options;
        map<string, IntegerOption> integer_options;
        // map<string, StringOption> string_options;
        map<string, BooleanOption> boolean_options;
    };

} // namespace fatrop
#endif // FatropParams