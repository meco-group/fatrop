/*
 * Fatrop - A fast trajectory optimization solver
 *  Copyright (C) 2022 - 2024 Lander Vanroye, KU Leuven. All rights reserved.
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
#ifndef __fatrop_solver_FatropOptions_hpp__
#define __fatrop_solver_FatropOptions_hpp__
#include <iostream>
#include <string>
#include <limits>
#include <unordered_map>
#include <memory>
#include <cstring>

namespace fatrop
{
    /**
     * @enum OptionType
     *
     * Enum that defines the different types of options supported by the system.
     */
    enum OptionType
    {
        INT,    ///< Represents an integer option.
        DOUBLE, ///< Represents a double-precision floating-point option.
        BOOL,   ///< Represents a boolean option.
        STRING  ///< Represents a string option.
    };
    // forward declaration
    struct OptionValueVariant;

    /**
     * @class OptionBase
     *
     * Base class for different option types. Provides a generic interface
     * for managing option values of various types.
     */
    class OptionBase
    {
    public:
        OptionBase(const std::string &name, const std::string &description) : name(name), description(description) {}
        virtual void set(int value) { throw std::invalid_argument("Invalid type for option " + name); };
        virtual void set(double value) { throw std::invalid_argument("Invalid type for option " + name); };
        virtual void set(const std::string &value) { throw std::invalid_argument("Invalid type " + name); };
        virtual void set(bool value) { throw std::invalid_argument("Invalid type " + name); };
        virtual void set(const OptionValueVariant &value);
        std::string name;
        std::string description;
    };

    /**
     * @class Option<T>
     *
     * Templated class that holds a value of type T and provides methods
     * to set and get that value. Inherits from OptionBase.
     *
     * @tparam T The type of the option (e.g., int, double, std::string).
     */
    template <typename T>
    class Option : public OptionBase
    {
    public:
        /**
         * @brief Constructor that initializes the option with a name, description, and default value.
         * @param name The name of the option.
         * @param description The description of the option.
         * @param default_value The default value of the option.
         */
        Option(const ::std::string &name, const std::string &description, const T &default_value) : OptionBase(name, description), value(default_value) {}
        void set_value(T value_in) { this->value = value_in; }
        const T &get() const { return value; }

    private:
        T value;
    };
    /**
     * @class BoolOption
     *
     * A subclass of Option<bool> for handling boolean options. Allows setting
     * values using either a boolean or a string ("yes"/"no").
     */
    class BoolOption : public Option<bool>
    {
    public:
        BoolOption(const std::string &name, const std::string &description, bool default_value) : Option<bool>(name, description, default_value) {}
        void set(bool value) { Option<bool>::set_value(value); }
        void set(const std::string &value);
    };

    /**
     * @class NumericOption<T>
     *
     * A subclass of Option<T> that handles numeric options, such as int or double,
     * with the ability to enforce lower and/or upper bounds.
     *
     * @tparam T The numeric type of the option (e.g., int, double).
     */
    template <typename T>
    class NumericOption : public Option<T>
    {
    public:
        /**
         * @brief Creates a lower-bounded numeric option.
         * @param name The name of the option.
         * @param description The description of the option.
         * @param default_value The default value of the option.
         * @param lower_bound The lower bound for the option value.
         * @return A new NumericOption with the lower bound.
         */
        static NumericOption<T> lower_bounded(const std::string &name, const std::string &description, const T &default_value, const T &lower_bound);
        /**
         * @brief Creates an upper-bounded numeric option.
         * @param name The name of the option.
         * @param description The description of the option.
         * @param default_value The default value of the option.
         * @param upper_bound The upper bound for the option value.
         * @return A new NumericOption with the upper bound.
         */
        static NumericOption<T> upper_bounded(const std::string &name, const std::string &description, const T &default_value, const T &upper_bound);
        /**
         * @brief Creates a bounded numeric option with both lower and upper bounds.
         * @param name The name of the option.
         * @param description The description of the option.
         * @param default_value The default value of the option.
         * @param lower_bound The lower bound for the option value.
         * @param upper_bound The upper bound for the option value.
         * @return A new NumericOption with both bounds.
         */
        static NumericOption<T> bounded(const std::string &name, const std::string &description, const T &default_value, const T &lower_bound, const T &upper_bound);
        void set(T value_in);

    private:
        /**
         * @brief Constructor that initializes a numeric option with bounds.
         * @param name The name of the option.
         * @param description The description of the option.
         * @param default_value The default value of the option.
         * @param lower_bounded Whether the option has a lower bound.
         * @param lower_bound The lower bound.
         * @param upper_bounded Whether the option has an upper bound.
         * @param upper_bound The upper bound.
         */
        NumericOption(const std::string &name, const std::string &description, const T &default_value, bool lower_bounded, const T &lower_bound, bool upper_bounded, const T &upper_bound) : Option<T>(name, description, default_value), is_lower_bounded(lower_bounded), lower_bound(lower_bound), is_upper_bounded(upper_bounded), upper_bound(upper_bound) {}
        bool is_lower_bounded; ///< Whether the option has a lower bound.
        T lower_bound;         ///< The lower bound for the option value.
        bool is_upper_bounded; ///< Whether the option has an upper bound.
        T upper_bound;         ///< The upper bound for the option value.
    };

    typedef NumericOption<int> IntOption;
    typedef NumericOption<double> DoubleOption;
    typedef Option<std::string> StringOption;
    typedef BoolOption BoolOption;

    /**
     * @struct FatropOptions
     *
     * Contains the configuration settings for the system.
     * It is optimized for fast look-up performance,
     * as these options are accessed each time the Fatrop algorithm is executed.
     */
    struct FatropOptions
    {
        // print options
        IntOption print_level = IntOption::lower_bounded("print_level", "fatrop print level", 10, 0);
        // fatrop algorithm options
        IntOption max_iter = IntOption::lower_bounded("max_iter", "maximum number of iterations", 1000, 0);
        DoubleOption tol = DoubleOption::lower_bounded("tol", "tolerance", 1e-8, 0.0);
        DoubleOption acceptable_tol = DoubleOption::lower_bounded("acceptable_tol", "acceptable tolerance", 1e-6, 0.0);
        IntOption max_watchdog_steps = IntOption::lower_bounded("max_watchdog_steps", "maximum number of watchdog steps", 4, 0);
        IntOption acceptable_iter = IntOption::lower_bounded("acceptable_iter", "acceptable iter", 15, 0);
        DoubleOption lammax = DoubleOption::lower_bounded("lammax", "lammax", 1e3, 0.0);
        DoubleOption mu_init = DoubleOption::lower_bounded("mu_init", "mu_init", 1e2, 0.0);
        DoubleOption kappa_eta = DoubleOption::lower_bounded("kappa_eta", "kappa_eta", 10.0, 0.0);
        DoubleOption kappa_mu = DoubleOption::lower_bounded("kappa_mu", "kappa_mu", 0.2, 0.0);
        DoubleOption theta_mu = DoubleOption::lower_bounded("theta_mu", "theta_mu", 1.5, 0.0);
        DoubleOption delta_w0 = DoubleOption::lower_bounded("delta_w0", "delta_w0", 1e-4, 0.0);
        DoubleOption delta_wmin = DoubleOption::lower_bounded("delta_wmin", "delta_wmin", 1e-20, 0.0);
        DoubleOption kappa_wmin = DoubleOption::lower_bounded("kappa_wmin", "kappa_wmin", 1.0 / 3.0, 0.0);
        DoubleOption kappa_wplus = DoubleOption::lower_bounded("kappa_wplus", "kappa_wplus", 8.0, 0.0);
        DoubleOption kappa_wplusem = DoubleOption::lower_bounded("kappa_wplusem", "kappa_wplusem", 100.0, 0.0);
        DoubleOption delta_c_stripe = DoubleOption::lower_bounded("delta_c_stripe", "delta_c_stripe", 1e-6, 0.0);
        DoubleOption kappa_c = DoubleOption::lower_bounded("kappa_c", "kappa_c", 0.25, 0.0);
        BoolOption warm_start_init_point = BoolOption("warm_start_init_point", "warm_start_init_point", false);
        DoubleOption theta_min = DoubleOption::lower_bounded("theta_min", "theta_min", 1e-4, 0.0);
        BoolOption recalc_y = BoolOption("recalc_y", "recalc_y", false);
        DoubleOption recalc_y_feas_tol = DoubleOption::lower_bounded("recalc_y_feas_tol", "recalc_y_feas_tol", 1e-6, 0.0);
        StringOption inequality_handling = StringOption("inequality_handling", "inequality_handling", "pd_ip");
        DoubleOption kappa_d = DoubleOption::lower_bounded("kappa_d", "kappa_d", 1e-5, 0.0);
        // line search options
        BoolOption accept_every_trial_step = BoolOption("accept_every_trial_step", "accept every trial step", false);
        DoubleOption s_phi = DoubleOption::lower_bounded("s_phi", "s_phi", 2.3, 0.0);
        DoubleOption delta = DoubleOption::lower_bounded("delta", "delta", 1.0, 0.0);
        DoubleOption s_theta = DoubleOption::lower_bounded("s_theta", "s_theta", 1.1, 0.0);
        DoubleOption gamma_theta = DoubleOption::lower_bounded("gamma_theta", "gamma_theta", 1e-5, 0.0);
        DoubleOption gamma_phi = DoubleOption::lower_bounded("gamma_phi", "gamma_phi", 1e-8, 0.0);
        DoubleOption eta_phi = DoubleOption::lower_bounded("eta_phi", "eta_phi", 1e-8, 0.0);
        DoubleOption gamma_alpha = DoubleOption::lower_bounded("gamma_alpha", "gamma_alpha", 0.05, 0.0);
        IntOption max_soc = IntOption::lower_bounded("max_soc", "max_soc", 2, 0);
        // linear solver options
        BoolOption linsol_iterative_refinement = BoolOption("linsol_iterative_refinement", "iterative ref", true);
        BoolOption linsol_perturbed_mode = BoolOption("linsol_perturbed_mode", "linear solver perturbed mode", false);
        BoolOption linsol_diagnostic = BoolOption("linsol_diagnostic", "linear solver diagnostic mode", false);
        DoubleOption linsol_perturbed_mode_param = DoubleOption::lower_bounded("linsol_perturbed_mode_param", "linear solver perturbed mode param", 1e-6, 0.);
        IntOption linsol_min_it_ref = IntOption::lower_bounded("linsol_min_it_ref", "minimum number of iterative refinement steps", 0, 0);
        IntOption linsol_max_it_ref = IntOption::lower_bounded("linsol_max_it_ref", "maximum number of iterative refinement steps", 5, 0);
        DoubleOption linsol_min_it_acc = DoubleOption::lower_bounded("linsol_min_it_acc", "stopping criterion for iterative refinement procedure", 1e-8, 0.);
        DoubleOption linsol_lu_fact_tol = DoubleOption::lower_bounded("linsol_lu_fact_tol", "pivoting tolerance parameter for lu fact", 1e-5, 0.);
        BoolOption iterative_refinement_SOC = BoolOption("iterative_refinement_SOC", "Use iterative refinement for SOC", true);
        BoolOption ls_scaling = BoolOption("ls_scaling", "Use automatic scaling for linear system", true);
        // fatrop data options
        DoubleOption warm_start_mult_bound_push = DoubleOption::lower_bounded("warm_start_mult_bound_push", "warm_start_mult_bound_push", 1e-2, 0.0);
        DoubleOption smax = DoubleOption::lower_bounded("smax", "smax", 100.0, 0.0);
        DoubleOption bound_push = DoubleOption::lower_bounded("bound_push", "kappa1", 1e-2, 0.0);
        DoubleOption bound_frac = DoubleOption::lower_bounded("bound_frac", "kappa2", 1e-2, 0.0);
        DoubleOption kappa_sigma = DoubleOption::lower_bounded("kappa_sigma", "kappa_sigma", 1e10, 0.0);
        DoubleOption bound_relax_factor = DoubleOption::lower_bounded("bound_relax_factor", "bound_relax_factor", 1e-8, 0.0);
        DoubleOption constr_viol_tol = DoubleOption::lower_bounded("constr_viol_tol", "constr_viol_tol", 1e-4, 0.0);
        // restoration phase options
        DoubleOption resto_rho = DoubleOption::lower_bounded("resto_rho", "Resto L1 penalty parameter", 1000., 0.0);
        DoubleOption resto_xi = DoubleOption::lower_bounded("resto_xi", "Resto xi parameter", 1., 0.0);
    };

    /**
     * @struct OptionValueVariant
     *
     * A structure used to hold a value of different types (int, double, string)
     * and identify its type. This structure is used to pass values to options.
     */
    struct OptionValueVariant
    {
        OptionValueVariant(int value_in) : type(OptionType::INT), value(malloc(sizeof(int))) { *static_cast<int *>(value) = value_in; }
        OptionValueVariant(double value_in) : type(OptionType::DOUBLE), value(malloc(sizeof(double))) { *static_cast<double *>(value) = value_in; }
        OptionValueVariant(bool value_in) : type(OptionType::BOOL), value(malloc(sizeof(bool))) { *static_cast<bool *>(value) = value_in; }
        OptionValueVariant(const std::string &value_in) : type(OptionType::STRING), value(malloc((value_in.size()+1)*sizeof(char))) { std::strcpy(static_cast<char *>(value), value_in.c_str()); }
        OptionValueVariant(const char *value_in) : type(OptionType::STRING), value(malloc((std::strlen(value_in)+1)*sizeof(char))) { std::strcpy(static_cast<char *>(value), value_in); }
        ~OptionValueVariant() { free(value); }
        OptionType type;
        void *value;
    };

    /**
     * @class FatropOptionsRegistry
     *
     * A class that manages a collection of options (FatropOptions).
     * Provides a way to set option values dynamically.
     */
    class FatropOptionsRegistry
    {
    public:
        FatropOptionsRegistry(FatropOptions &options);
        // T can be int, double, bool, string, or an OptionValueVariant
        template <typename T>
        void set(const std::string &name, T value);
    public:
        std::unordered_map<std::string, OptionBase *> options;
    };

};

#endif // __fatrop_solver_FatropOptions_hpp__