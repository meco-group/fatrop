//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//
#ifndef __fatrop_ip_iteration_output_hxx__
#define __fatrop_ip_iteration_output_hxx__
#include "fatrop/common/printing.hpp"
#include "fatrop/ip_algorithm/ip_iteration_output.hpp"
#include <cmath>
#include <iomanip>
#include <string>

namespace fatrop
{
    template <typename ProblemType>
    IpIterationOutput<ProblemType>::IpIterationOutput(const IpDataSp &ipdata) : ipdata_(ipdata)
    {
    }

    template <typename ProblemType> void IpIterationOutput<ProblemType>::print_header()
    {
        f_out << std::setw(4) << "iter" << " " << std::setw(12) << "objective" << " "
                   << std::setw(8) << "inf_pr" << " " << std::setw(8) << "inf_du" << " "
                   << std::setw(6) << "lg(mu)" << " " << std::setw(8) << "||d||" << " "
                   << std::setw(8) << "lg(rg)"
                   << " " << std::setw(10) << "alpha_du" << " " << std::setw(10) << "alpha_pr"
                   << " " << std::setw(4) << "ls" << std::endl;
    }

    template <typename ProblemType>
    void IpIterationOutput<ProblemType>::print_iteration(Index iter, Scalar objective,
                                                         Scalar inf_pr, Scalar inf_du, Scalar lg_mu,
                                                         Scalar d_norm, Scalar rg, Scalar alpha_du,
                                                         Scalar alpha_pr, Index ls)
    {
        f_out << std::setw(4) << iter << " " << std::setw(12) << std::scientific
                   << std::setprecision(8) << objective << " " << std::setw(8) << std::scientific
                   << std::setprecision(2) << inf_pr << " " << std::setw(8) << std::scientific
                   << std::setprecision(2) << inf_du << " " << std::setw(6) << std::fixed
                   << std::setprecision(1) << lg_mu << " " << std::setw(8) << std::scientific
                   << std::setprecision(2) << d_norm << " " << std::setw(6);
        if (rg == 0.0)
            f_out << "-";
        else
            f_out << std::fixed << std::setprecision(1) << std::log10(rg);
        f_out << std::setw(10) << std::scientific << std::setprecision(2) << alpha_du << " "
                   << std::setw(10) << std::scientific << std::setprecision(2) << alpha_pr << " "
                   << std::setw(4) << ls << std::endl;
    }

    template <typename ProblemType> void IpIterationOutput<ProblemType>::output_current_iteration()
    {
        IpIterateType &curr_it = ipdata_->current_iterate();
        const Index iter = ipdata_->iteration_number();
        const Scalar objective = curr_it.obj_value();
        const Scalar inf_pr = norm_inf(curr_it.constr_viol());
        const Scalar inf_du = std::max(norm_inf(curr_it.dual_infeasibility_x()),
                                       norm_inf(curr_it.dual_infeasibility_s()));
        const Scalar mu = curr_it.mu();
        const Scalar d_norm = curr_it.step_info().step_length;
        const Scalar rg = curr_it.step_info().inertia_correction_primal;
        const Scalar alpha_pr = curr_it.step_info().alpha_primal;
        const Scalar alpha_du = curr_it.step_info().alpha_dual;
        const Scalar ls = curr_it.step_info().ls_iter + 1;
        print_iteration(iter, objective, inf_pr, inf_du, std::log10(mu), d_norm, rg, alpha_du,
                        alpha_pr, ls);
    }
} // namespace fatrop

#endif // __fatrop_ip_iteration_output_hxx__
