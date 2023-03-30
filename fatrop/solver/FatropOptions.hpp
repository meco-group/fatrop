#ifndef FATROPPARAMSINCLUDED
#define FATROPPARAMSINCLUDED
#include "aux/SmartPtr.hpp"
#include <string>
#include <map>
#include <string>
#include <type_traits>
#include <iostream>
using namespace std;
namespace fatrop
{
    template <typename T>
    struct Option
    {
    public:
        // NumericOption operator=(const NumericOption &other) = default;
        Option(){};
        Option(const string &name, const string &description, T *value, T default_value, bool lower_bound_inclusive, T lower_bound, bool upper_bound_inclusive, T upper_bound); 
        static Option<T> LowerBounded(const string &name, const string &description, T *value, T default_value, T lower_bound);
        static Option<T> UpperBounded(const string &name, const string &description, T *value, T default_value, T upper_bound);
        static Option<T> UnBounded(const string &name, const string &description, T *value, T default_value);
        static Option<T> BoxBounded(const string &name, const string &description, T *value, T default_value, T lower_bound, T upper_bound);
        void SetDefault() const;
        void set(const T &new_value);
        string name_;
        string description_;
        T *value = NULL;
        T default_value_;
        bool lower_bound_inclusive_;
        T lower_bound_;
        bool upper_bound_inclusive_;
        T upper_bound_;
    };
    template <>
    struct Option<bool>
    {
        Option(){};
        Option(const string &name, const string &description, bool *value, bool default_value);
        void set(const bool &new_value);
        void SetDefault() const;
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
        FatropOptions();
        // the following options are shared between different algorithm components:
        int maxiter = 1000; // TODO this value cannot be changed to a value larger than the one used for building the solver
        double kappa_d = 1e-5;
        template <typename T>
        void set(const string &option_name, T value);

    public:
        void register_option(const NumericOption &option);
        void register_option(const IntegerOption &option);
        void register_option(const BooleanOption &option);
        friend auto operator<<(std::ostream &os, const FatropOptions &m) -> std::ostream &;
        map<string, NumericOption> numeric_options;
        map<string, IntegerOption> integer_options;
        // map<string, StringOption> string_options;
        map<string, BooleanOption> boolean_options;
    };

} // namespace fatrop
#endif // FatropOptions