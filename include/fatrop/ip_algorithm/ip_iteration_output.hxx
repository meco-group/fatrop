//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//
#ifndef __fatrop_ip_iteration_output_hxx__
#define __fatrop_ip_iteration_output_hxx__
#include "fatrop/ip_algorithm/ip_iteration_output.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace fatrop
{
    template <typename ProblemType>
    IpIterationOutput<ProblemType>::IpIterationOutput(const IpDataSp &ipdata) : ipdata_(ipdata)
    {
    }

    template <typename ProblemType> void IpIterationOutput<ProblemType>::print_header()
    {
        std::cout << std::setw(4) << "iter" << " " << std::setw(12) << "objective" << " "
                  << std::setw(8) << "inf_pr" << " " << std::setw(8) << "inf_du" << " "
                  << std::setw(6) << "lg(mu)" << " " << std::setw(8) << "||d||" << " "
                  << std::setw(8) << "lg(rg)"
                  << " " << std::setw(10) << "alpha_du" << " " << std::setw(10) << "alpha_pr" << " "
                  << std::setw(4) << "ls" << std::endl;
    }

    template <typename ProblemType>
    void IpIterationOutput<ProblemType>::print_iteration(Index iter, Scalar objective,
                                                         Scalar inf_pr, Scalar inf_du, Scalar lg_mu,
                                                         Scalar d_norm, Scalar rg, Scalar alpha_du,
                                                         Scalar alpha_pr, Index ls)
    {
        std::string lg_rg = (rg == 0.0) ? "-" : std::to_string(std::log10(rg));

        std::cout << std::setw(4) << iter << " " << std::setw(12) << std::scientific
                  << std::setprecision(8) << objective << " " << std::setw(8) << std::scientific
                  << std::setprecision(2) << inf_pr << " " << std::setw(8) << std::scientific
                  << std::setprecision(2) << inf_du << " " << std::setw(6) << std::fixed
                  << std::setprecision(1) << lg_mu << " " << std::setw(8) << std::scientific
                  << std::setprecision(2) << d_norm << " " << std::setw(6) << lg_rg << " "
                  << std::setw(10) << std::scientific << std::setprecision(2) << alpha_du << " "
                  << std::setw(10) << std::scientific << std::setprecision(2) << alpha_pr << " "
                  << std::setw(4) << ls << std::endl;
    }

    template <typename ProblemType> void IpIterationOutput<ProblemType>::output_current_iteration()
    {
        IpIterateType &curr_it = ipdata_->current_iterate();
        const Index iter = ipdata_->iteration_number();
        const Scalar objective = curr_it.obj_value();
        const Scalar inf_pr = norm_inf(curr_it.constr_viol());
        const Scalar inf_du = std::max(norm_inf(curr_it.dual_infeasibility_x()), norm_inf(curr_it.dual_infeasibility_s()));
        const Scalar mu = curr_it.mu();
        const Scalar d_norm = std::max(norm_inf(curr_it.delta_primal_x()), norm_inf(curr_it.delta_primal_s()));
        const Scalar rg = 0.;
        const Scalar alpha_pr = 0.;
        const Scalar alpha_du = 0.;
        const Scalar ls = 0.;
        print_iteration(iter, objective, inf_pr, inf_du, std::log10(mu), d_norm, rg, alpha_du, alpha_pr, ls);
    }
} // namespace fatrop

#endif // __fatrop_ip_iteration_output_hxx__
