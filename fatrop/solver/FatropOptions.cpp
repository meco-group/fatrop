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
#include "FatropOptions.hpp"
using namespace fatrop;

// void OptionBase::set(const OptionValueVariant &value)
// {
//     switch (value.type)
//     {
//     case OptionType::INT:
//         set(*static_cast<int *>(value.value));
//         break;
//     case OptionType::DOUBLE:
//         set(*static_cast<double *>(value.value));
//         break;
//     case OptionType::BOOL:
//         set(*static_cast<bool *>(value.value));
//         break;
//     case OptionType::STRING:
//         set(std::string(static_cast<char *>(value.value)));
//         break;
//     default:
//         throw std::invalid_argument("Invalid type for option " + name);
//     }
// }
void BoolOption::set(const std::string &value)
{
    if (value == "yes")
    {
        Option<bool>::set_value(true);
    }
    else if (value == "no")
    {
        Option<bool>::set_value(false);
    }
    else
    {
        throw std::invalid_argument("Invalid string value for boolean option, should be 'yes' or 'no'");
    }
}
template <typename T>
NumericOption<T> NumericOption<T>::lower_bounded(const std::string &name, const std::string &description, const T &default_value, const T &lower_bound, bool requires_reinitialization)
{
    return NumericOption<T>(name, description, default_value, true, lower_bound, false, 0, requires_reinitialization);
}
template <typename T>
NumericOption<T> NumericOption<T>::upper_bounded(const std::string &name, const std::string &description, const T &default_value, const T &upper_bound, bool requires_reinitialization)
{
    return NumericOption<T>(name, description, default_value, false, 0, true, upper_bound, requires_reinitialization);
}
template <typename T>
NumericOption<T> NumericOption<T>::bounded(const std::string &name, const std::string &description, const T &default_value, const T &lower_bound, const T &upper_bound, bool requires_reinitialization)
{
    return NumericOption<T>(name, description, default_value, true, lower_bound, true, upper_bound, requires_reinitialization);
}
template <typename T>
void NumericOption<T>::set(T value_in)
{
    if ((is_lower_bounded && (value_in < lower_bound)) || (is_upper_bounded && (value_in > upper_bound)))
    {
        throw std::invalid_argument("Value out of bounds");
    }
    Option<T>::set_value(value_in);
}
// instantiate the template for int and double
template class NumericOption<int>;
template class NumericOption<double>;

FatropOptionsRegistry::FatropOptionsRegistry(FatropOptions &options)
{
    // register all options
    OptionBase *all_options[] = {
        &options.max_iter,
        &options.print_level,
        &options.tol,
        &options.acceptable_tol,
        &options.max_watchdog_steps,
        &options.acceptable_iter,
        &options.lammax,
        &options.mu_init,
        &options.kappa_eta,
        &options.kappa_mu,
        &options.theta_mu,
        &options.delta_w0,
        &options.delta_wmin,
        &options.kappa_wmin,
        &options.kappa_wplus,
        &options.kappa_wplusem,
        &options.delta_c_stripe,
        &options.kappa_c,
        &options.warm_start_init_point,
        &options.theta_min,
        &options.recalc_y,
        &options.recalc_y_feas_tol,
        &options.inequality_handling,
        &options.kappa_d,
        &options.accept_every_trial_step,
        &options.s_phi,
        &options.delta,
        &options.s_theta,
        &options.gamma_theta,
        &options.gamma_phi,
        &options.eta_phi,
        &options.gamma_alpha,
        &options.max_soc,
        &options.linsol_iterative_refinement,
        &options.linsol_perturbed_mode,
        &options.linsol_diagnostic,
        &options.linsol_perturbed_mode_param,
        &options.linsol_min_it_ref,
        &options.linsol_max_it_ref,
        &options.linsol_min_it_acc,
        &options.linsol_lu_fact_tol,
        &options.iterative_refinement_SOC,
        &options.ls_scaling,
        &options.warm_start_mult_bound_push,
        &options.smax,
        &options.bound_push,
        &options.bound_frac,
        &options.kappa_sigma,
        &options.bound_relax_factor,
        &options.constr_viol_tol};
    for (auto option : all_options)
    {
        this->options[option->name] = option;
    }
}

template<typename T>
void FatropOptionsRegistry::set(const std::string &name, T value)
{
    if (options.find(name) == options.end())
    {
        throw std::invalid_argument("Option with name " + name + " not found");
    }
    options[name]->set(value);
}

template void FatropOptionsRegistry::set<int>(const std::string &, int);
template void FatropOptionsRegistry::set<double>(const std::string &, double);
template void FatropOptionsRegistry::set<bool>(const std::string &, bool);
template void FatropOptionsRegistry::set<std::string>(const std::string &, std::string);
template void FatropOptionsRegistry::set<const std::string&>(const std::string &, const std::string&);
// template void FatropOptionsRegistry::set<OptionValueVariant>(const std::string &, OptionValueVariant);
// template void FatropOptionsRegistry::set<const OptionValueVariant&>(const std::string &, const OptionValueVariant&);