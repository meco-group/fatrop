#ifndef FATROPPARAMSINCLUDED
#define FATROPPARAMSINCLUDED
#include "aux/SmartPtr.hpp"
#include <string>
#include <map>
#include <string>
using namespace std;
namespace fatrop
{
    struct NumericOption
    {
    public:
        // NumericOption operator=(const NumericOption &other) = default;
        NumericOption(){};
        NumericOption(const string &name, const string &description, double *value, double default_value, bool lower_bound_inclusive, double lower_bound, bool upper_bound_inclusive, double upper_bound):
            name_(name), description_(description), value(value), default_value_(default_value), lower_bound_inclusive_(lower_bound_inclusive), lower_bound_(lower_bound), upper_bound_inclusive_(upper_bound_inclusive), upper_bound_(upper_bound)
        {
        };
        static NumericOption LowerBounded(const string &name, const string &description, double *value, double default_value, double lower_bound)
        {
            return NumericOption(name, description, value, default_value, true, lower_bound, false, 0.0);
        };
        static NumericOption UpperBounded(const string &name, const string &description, double *value, double default_value, double upper_bound)
        {
            return NumericOption(name, description, value, default_value, false, 0.0, true, upper_bound);
        };
        static NumericOption UnBounded(const string &name, const string &description, double *value, double default_value)
        {
            return NumericOption(name, description, value, default_value, false, 0.0, false, 0.0);
        };
        static NumericOption BoxBounded(const string &name, const string &description, double *value, double default_value, double lower_bound, double upper_bound)
        {
            return NumericOption(name, description, value, default_value, true, lower_bound, true, upper_bound);
        };
        string name_;
        string description_;
        double *value = NULL;
        double default_value_;
        bool lower_bound_inclusive_;
        double lower_bound_;
        bool upper_bound_inclusive_;
        double upper_bound_;
    };

    class FatropOptions
    {
    public:
        FatropOptions()
        {
            // register tolerance option
            RegisterNumericOption(NumericOption::LowerBounded("tol", "tolerance", &tol, 1e-8, 0.0));
            RegisterNumericOption(NumericOption::LowerBounded("acceptable_tol", "acceptable tolerance", &acceptable_tol, 1e-6, 0.0));
        };

        int max_watchdog_steps = 4;
        bool warm_start_dual = false;
        int maxiter = 1000; // TODO unsafe to change maxiter because it it used at building!!
        double smax = 100.0;
        double lammax = 1e3;
        double tol = 1e-8;
        double acceptable_tol = 1e-6;
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
        double kappa_d = 1e-5;
        double bound_relax_factor = 1e-8;
        double constr_viol_tol = 1e-4; // currently only used to relax bounds
        void SetNumericOption(const string& option_name, double value)
        {
            // check if option exists
            if (numeric_options.find(option_name) == numeric_options.end())
            {
                throw runtime_error("Option " + option_name + " does not exist.");
            }
            // check if value is in bounds
            NumericOption& option = numeric_options[option_name];
            if (option.lower_bound_inclusive_ && value < option.lower_bound_)
            {
                throw runtime_error("Option " + option_name + " must be greater than or equal to " + to_string(option.lower_bound_));
            }
            if (option.upper_bound_inclusive_ && value > option.upper_bound_)
            {
                throw runtime_error("Option " + option_name + " must be less than or equal to " + to_string(option.upper_bound_));
            }
            // set value
            *option.value = value;
        }
    private:
        void RegisterNumericOption(const NumericOption &option)
        {
            numeric_options[option.name_] = option;
        }
        map<string, NumericOption> numeric_options;
    };

} // namespace fatrop
#endif // FatropParams