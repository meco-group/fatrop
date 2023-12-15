#include "ustage_quantities.hpp"
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
#include "fatrop/spectool/auxiliary/constraint_helper.hpp"
#include "ustage.hpp"
#include "ocp.hpp"
namespace fatrop
{
    namespace spectool
    {
        uStageQuantities uStageQuantities::create(const std::shared_ptr<const uStageInternal> &ustage, const std::shared_ptr<const uStageInternal> &prev, const std::shared_ptr<const uStageInternal> &next, const std::vector<cs::MX> &global_parameter_syms)
        {
            uStageQuantities ret;
            ret.x = cs::MX::veccat(ustage->get_states(true, prev));
            ret.p_stage = cs::MX::veccat(ustage->get_control_parameters());
            ret.u = cs::MX::veccat(ustage->get_controls(true, prev));
            ret.p_global = cs::MX::veccat(global_parameter_syms);
            ret.K = ustage->K();
            // std::cout << "number of hybrids that are states " << ustage->get_hybrids_states(prev).size() << std::endl;
            // std::cout << "number of hybrids that are controls " << ustage->get_hybrids_controls(prev).size() << std::endl;
            cs::MX x_next;
            if (next)
            {
                std::vector<cs::MX> from;
                std::vector<cs::MX> to;
                for (auto &[from_, to_] : ustage->get_next_states())
                {
                    from.push_back(from_);
                    to.push_back(to_);
                }
                // check if every state of the next stage is defined
                for (auto &state : next->get_states(false))
                {
                    if (ustage->get_next_states().find(state) == ustage->get_next_states().end())
                    {
                        throw std::runtime_error("Did you set_next for every state?");
                    }
                }
                try
                {
                    auto x_next_vec = next->get_states(true, ustage);
                    ret.xp1 = cs::MX::sym("xp1", cs::MX::veccat(next->get_states(true, ustage)).size1());
                    x_next = cs::MX::veccat(cs::MX::substitute({x_next_vec}, from, to));
                }
                catch (const std::exception &e)
                {
                    std::cerr << e.what() << '\n';
                    throw std::runtime_error("Did you set_next for every state?");
                }
                if (ustage->K() > 1)
                {
                    // check if x_next is same as x
                    auto x_next_vec = next->get_states(true, ustage);
                    auto x_vec = ustage->get_states(true, prev);
                    // check if size is same
                    if (x_next_vec.size() != x_vec.size())
                    {
                        throw std::runtime_error("x_next and x must have the same size");
                    }
                    // iterate over elements and check if the same
                    for (int i = 0; i < x_next_vec.size(); i++)
                    {
                        if (x_next_vec[i].get() != x_vec[i].get())
                        {
                            throw std::runtime_error("x_next and x must be the same");
                        }
                    }
                }
                // ret.nxp1 = ret.x_next.size1();
            }
            else
            {
                // ret.nxp1 = 0;
                x_next = cs::MX::zeros(0, 1);
                ret.xp1 = cs::MX::zeros(0, 1);
            }
            cs::MX g;
            cs::MX g_ineq;
            ConstraintHelper::process(ustage->get_constraints(),
                                      ret.lb,
                                      ret.ub,
                                      g_ineq,
                                      g);
            ret.lam_dyn = cs::MX::sym("lam_dyn", x_next.size1(), 1);
            ret.lam_g_ineq = cs::MX::sym("lam_g_ineq", g_ineq.size1(), 1);
            ret.lam_g_eq = cs::MX::sym("lam_g_eq", g.size1(), 1);
            auto ux = cs::MX::veccat({ret.u, ret.x});
            ret.Gg_dyn = uStageQuantities::dynamics_jacobian_sym(ux, ret.xp1, x_next);
            ret.Gg_ineq = uStageQuantities::inequality_jacobian_sym(ux, g_ineq);
            ret.Gg_eq = uStageQuantities::equality_jacobian_sym(ux, g);



            ret.L = 0;
            for (auto &term : ustage->get_objective_terms())
            {
                ret.L += term;
            }
            ret.hess_obj = uStageQuantities::hess_lag_obj_sym(ux, ret.L);
            ret.hess_dyn = uStageQuantities::hess_lag_dyn_sym(ux, x_next, ret.lam_dyn);
            ret.hess_g_ineq = uStageQuantities::hess_lag_dyn_sym(ux, g_ineq, ret.lam_g_ineq);
            ret.hess_g_eq = uStageQuantities::hess_lag_dyn_sym(ux, g, ret.lam_g_eq);
            return ret;
        }
    }
}