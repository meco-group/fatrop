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
        void SetDefault()
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
            // register integer options
            RegisterOption(IntegerOption::BoxBounded("max_iter", "maximum number of iterations", &maxiter, 1000, 0, maxiter));
            RegisterOption(IntegerOption::LowerBounded("max_watchdog_steps", "maximum number of watchdog steps", &max_watchdog_steps, 4, 0));
            // RegisterIntegerOption(IntegerOption::LowerBounded("max_watchdog_steps", "maximum number of watchdog steps", &max_watchdog_steps, 4, 0));
            // register tolerance option
            RegisterOption(NumericOption::LowerBounded("tol", "tolerance", &tol, 1e-8, 0.0));
            RegisterOption(NumericOption::LowerBounded("acceptable_tol", "acceptable tolerance", &acceptable_tol, 1e-6, 0.0));
            RegisterOption(NumericOption::LowerBounded("smax", "smax", &smax, 100.0, 0.0));
            RegisterOption(NumericOption::LowerBounded("lammax", "lammax", &lammax, 1e3, 0.0));
            RegisterOption(NumericOption::LowerBounded("mu_init", "mu_init", &mu0, 1e2, 0.0));
            RegisterOption(NumericOption::LowerBounded("kappa_eta", "kappa_eta", &kappa_eta, 10.0, 0.0));
            RegisterOption(NumericOption::LowerBounded("kappa_mu", "kappa_mu", &kappa_mu, 0.2, 0.0));
            RegisterOption(NumericOption::LowerBounded("theta_mu", "theta_mu", &theta_mu, 1.5, 0.0));
            RegisterOption(NumericOption::LowerBounded("delta_w0", "delta_w0", &delta_w0, 1e-4, 0.0));
            RegisterOption(NumericOption::LowerBounded("delta_wmin", "delta_wmin", &delta_wmin, 1e-20, 0.0));
            RegisterOption(NumericOption::LowerBounded("kappa_wmin", "kappa_wmin", &kappa_wmin, 1.0 / 3.0, 0.0));
            RegisterOption(NumericOption::LowerBounded("kappa_wplus", "kappa_wplus", &kappa_wplus, 8.0, 0.0));
            RegisterOption(NumericOption::LowerBounded("kappa_wplusem", "kappa_wplusem", &kappa_wplusem, 100.0, 0.0));
            RegisterOption(NumericOption::LowerBounded("kappa_sigma", "kappa_sigma", &kappa_sigma, 1e10, 0.0));
            RegisterOption(NumericOption::LowerBounded("s_phi", "s_phi", &s_phi, 2.3, 0.0));
            RegisterOption(NumericOption::LowerBounded("delta", "delta", &delta, 1.0, 0.0));
            RegisterOption(NumericOption::LowerBounded("s_theta", "s_theta", &s_theta, 1.1, 0.0));
            RegisterOption(NumericOption::LowerBounded("theta_min", "theta_min", &theta_min, 1e-4, 0.0));
            RegisterOption(NumericOption::LowerBounded("gamma_theta", "gamma_theta", &gamma_theta, 1e-12, 0.0));
            RegisterOption(NumericOption::LowerBounded("gamma_phi", "gamma_phi", &gamma_phi, 1e-8, 0.0));
            RegisterOption(NumericOption::LowerBounded("gamma_alpha", "gamma_alpha", &gamma_alpha, 0.05, 0.0));
            RegisterOption(NumericOption::LowerBounded("eta_phi", "eta_phi", &eta_phi, 1e-8, 0.0));
            RegisterOption(NumericOption::LowerBounded("delta_c_stripe", "delta_c_stripe", &delta_c_stripe, 1e-2, 0.0));
            RegisterOption(NumericOption::LowerBounded("kappa_c", "kappa_c", &kappa_c, 0.25, 0.0));
            RegisterOption(NumericOption::LowerBounded("bound_push", "kappa1", &kappa1, 1e-2, 0.0));
            RegisterOption(NumericOption::LowerBounded("bound_frac", "kappa2", &kappa2, 1e-2, 0.0));
            RegisterOption(NumericOption::LowerBounded("warm_start_mult_bound_push", "warm_start_mult_bound_push", &warm_start_mult_bound_push, 1e-2, 0.0));
            RegisterOption(NumericOption::LowerBounded("kappa_d", "kappa_d", &kappa_d, 1e-5, 0.0));
            RegisterOption(NumericOption::LowerBounded("bound_relax_factor", "bound_relax_factor", &bound_relax_factor, 1e-8, 0.0));
            RegisterOption(NumericOption::LowerBounded("constr_viol_tol", "constr_viol_tol", &constr_viol_tol, 1e-4, 0.0));
            RegisterOption(BooleanOption("warm_start_init_point", "warm_start_init_point", &warm_start_init_point, false));
        };

        int max_watchdog_steps = 4;
        bool warm_start_init_point = false;
        int maxiter = 1000; // TODO this value cannot be changed to a value larger than the one used for building the solver
        double tol = 1e-8;
        double acceptable_tol = 1e-6;
        double smax = 100.0;
        double lammax = 1e3;
        double acceptable_iter = 15;
        double mu0 = 1e2;
        // double mu0 = 2e5;
        // double mu0 = 1e;
        double kappa_eta = 10;
        double kappa_mu = 0.2;
        double theta_mu = 1.5;
        double delta_w0 = 1e-4;
        double delta_wmin = 1e-20;
        double kappa_wmin = 1.0 / 3.0;
        double kappa_wplus = 8;
        double kappa_wplusem = 100;
        double kappa_sigma = 1e10;
        double s_phi = 2.3;
        double delta = 1.0;
        double s_theta = 1.1;
        double theta_min = 1e-4;
        // double gamma_theta = 1e-8;
        double gamma_theta = 1e-12;
        // double gamma_theta = 1e-5; // todo check!!
        double gamma_phi = 1e-8;
        double gamma_alpha = 0.05;
        double eta_phi = 1e-8;
        // double delta_c_stripe = 1e-8;
        double delta_c_stripe = 1e-2;
        double kappa_c = 0.25;
        double kappa1 = 1e-2;
        double kappa2 = 1e-2;
        double warm_start_mult_bound_push = 1e-2;
        double kappa_d = 1e-5;
        double bound_relax_factor = 1e-8;
        double constr_viol_tol = 1e-4; // currently only used to relax bounds
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
        }
        void RegisterOption(const IntegerOption &option)
        {
            integer_options[option.name_] = option;
        }
        // void RegisterStringOption(const StringOption &option)
        // {
        //     string_options[option.name_] = option;
        // }
        void RegisterOption(const BooleanOption &option)
        {
            boolean_options[option.name_] = option;
        }
        map<string, NumericOption> numeric_options;
        map<string, IntegerOption> integer_options;
        // map<string, StringOption> string_options;
        map<string, BooleanOption> boolean_options;
    };

} // namespace fatrop
#endif // FatropParams