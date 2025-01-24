#include "ip_search_dir.hpp"

#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_iterate.hpp"
#include <cmath>

namespace fatrop
{
    template <typename ProblemType, typename LinearSystemType, typename LinearSolverDerived>
    IpSearchDirImpl<ProblemType, LinearSystemType, LinearSolverDerived>::IpSearchDirImpl(
        const IpDataSp &ipdata, const LinearSolverSp &linear_solver)
        : ipdata_(ipdata), linear_solver_(linear_solver),
          rhs_x_(ipdata->current_iterate().nlp()->nlp_dims().number_of_variables),
          rhs_s_(ipdata->current_iterate().nlp()->nlp_dims().number_of_ineq_constraints),
          rhs_g_(ipdata->current_iterate().nlp()->nlp_dims().number_of_eq_constraints),
          rhs_cl_(ipdata->current_iterate().nlp()->nlp_dims().number_of_ineq_constraints),
          rhs_cu_(ipdata->current_iterate().nlp()->nlp_dims().number_of_ineq_constraints),
          Dx_(ipdata->current_iterate().nlp()->nlp_dims().number_of_variables),
          Ds_(ipdata->current_iterate().nlp()->nlp_dims().number_of_ineq_constraints),
          Deq_(ipdata->current_iterate().nlp()->nlp_dims().number_of_eq_constraints)
    {
    }
    template <typename ProblemType, typename LinearSystemType, typename LinearSolverDerived>
    void IpSearchDirImpl<ProblemType, LinearSystemType, LinearSolverDerived>::reset()
    {
        delta_w_last_ = 0.;
    }
    template <typename ProblemType, typename LinearSystemType, typename LinearSolverDerived>
    LinsolReturnFlag
    IpSearchDirImpl<ProblemType, LinearSystemType, LinearSolverDerived>::compute_search_dir()
    {
        IpIterateType &curr_it = ipdata_->current_iterate();
        rhs_x_.block(rhs_x_.m(), 0) = curr_it.dual_infeasibility_x();
        rhs_s_.block(rhs_s_.m(), 0) = curr_it.dual_infeasibility_s();
        rhs_g_.block(rhs_g_.m(), 0) = curr_it.constr_viol();
        rhs_cl_.block(rhs_cl_.m(), 0) = curr_it.delta_lower() * curr_it.dual_bounds_l();
        rhs_cu_.block(rhs_cu_.m(), 0) = curr_it.delta_upper() * curr_it.dual_bounds_u();

        Dx_ = 0.;
        Deq_ = 0.;

        Scalar delta_w = 0.;
        Scalar delta_c = 0.;
        bool perturbed_eq = false;
        Scalar mu = curr_it.mu();

        bool solved = false;
        bool first_try_delta_w = true;
        LinsolReturnFlag ret;
        while (!solved)
        {
            bool update_delta_w = false;
            bool update_delta_c = false;
            Dx_ = delta_w;
            Deq_ = delta_c;
            LinearSystem<LinearSystemType> ls(
                curr_it.info(), curr_it.jacobian(), curr_it.hessian(), Dx_, perturbed_eq, Deq_,
                curr_it.delta_lower(), curr_it.delta_upper(), curr_it.dual_bounds_l(),
                curr_it.dual_bounds_u(), rhs_x_, rhs_s_, rhs_g_, rhs_cl_, rhs_cu_);
            ret = linear_solver_->solve_in_place(ls);
            switch (ret)
            {
            case (LinsolReturnFlag::SUCCESS):
                solved = true;
                break;
            case (LinsolReturnFlag::ITREF_MAX_ITER):
                solved = true;
                break;
            case (LinsolReturnFlag::ITREF_INCREASE):
                solved = true;
                break;
            case (LinsolReturnFlag::INDEFINITE):
                update_delta_w = true;
                solved = false;
                break;
            case (LinsolReturnFlag::NOFULL_RANK):
                update_delta_c = true;
                perturbed_eq = true;
                solved = false;
                break;
            }

            if (update_delta_w)
            {
                if (first_try_delta_w)
                {
                    delta_w = (delta_w_last_ == 0.)
                                  ? delta_w0_
                                  : std::max(delta_w_last_ * kappa_wmin_, delta_wmin_);
                    first_try_delta_w = false;
                }
                else
                {
                    delta_w =
                        (delta_w_last_ == 0.) ? kappa_wplusem_ * delta_w : kappa_wplus_ * delta_w;
                }
            }
            if (update_delta_c)
            {
                delta_c = delta_c_stripe_ * pow(mu, kappa_c_);
            }

            if (delta_w > 0)
                delta_w_last_ = delta_w;
        }

        curr_it.set_delta_primal_x(rhs_x_);
        curr_it.set_delta_primal_s(rhs_s_);
        curr_it.set_delta_dual_eq(rhs_g_);
        curr_it.set_delta_dual_bounds_l(rhs_cl_);
        curr_it.set_delta_dual_bounds_u(rhs_cu_);
        return ret;
    }
} // namespace fatrop