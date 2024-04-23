#include "ustage_quantities.hpp"
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
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
            uo_map_mx<Hessian> hessians;
            uo_map_mx<Jacobian> jacobians;
            std::vector<cs::MX> x_next;
            if (next || ustage->K() > 1)
            {
                auto next_states = (ustage->K() == 1) ? next->get_states(true, ustage) : ustage->get_states(true, prev);
                for (auto &state : next_states)
                {
                    try
                    {
                        x_next.push_back(ustage->get_next_states().at(state));
                    }
                    catch (const std::exception &e)
                    {
                        std::cerr << e.what() << '\n';
                        throw std::runtime_error("Did you set_next for every state?");
                    }
                    hessians[x_next.back()] = ustage->get_next_state_hessians().at(state);
                    jacobians[x_next.back()] = ustage->get_next_state_jacobians().at(state);
                }
                if (next && ustage->K() > 1)
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
                    for (size_t i = 0; i < x_next_vec.size(); i++)
                    {
                        if (x_next_vec[i].get() != x_vec[i].get())
                        {
                            throw std::runtime_error("x_next and x must be the same");
                        }
                    }
                }
                const int nxp1 = cs::MX::veccat(x_next).size1();
                ret.xp1 = cs::MX::sym("x_next", nxp1, 1);
            }
            else
            {
                // ret.nxp1 = 0;
                x_next = {cs::MX::zeros(0, 1)};
                ret.xp1 = cs::MX::zeros(0, 1);
            }
            // process constraints
            std::vector<cs::MX> equality_constraints;
            std::vector<cs::MX> inequality_constraints;
            std::vector<cs::DM> lb_vec;
            std::vector<cs::DM> ub_vec;
            {
                auto constraints = ustage->get_constraints();
                for (auto &constraint : constraints)
                {
                    cs::MX g;
                    cs::MX g_ineq;
                    cs::DM lb;
                    cs::DM ub;
                    ConstraintHelper::process(constraint,
                                              lb,
                                              ub,
                                              g_ineq,
                                              g);
                    bool is_ineq = g_ineq.size1() > 0;
                    bool is_eq = g.size1() > 0;
                    if (is_ineq && is_eq)
                    {
                        throw std::runtime_error("Constraint must be either equality or inequality, try splitting up when using subject_to(constraint)");
                    }
                    else if (is_ineq)
                    {
                        inequality_constraints.push_back(g_ineq);
                        lb_vec.push_back(lb);
                        ub_vec.push_back(ub);
                        jacobians[g_ineq] = ustage->get_constraint_jacobians().at(constraint);
                        hessians[g_ineq] = ustage->get_constraint_hessians().at(constraint);
                    }
                    else if (is_eq)
                    {
                        equality_constraints.push_back(g);
                        jacobians[g] = ustage->get_constraint_jacobians().at(constraint);
                        hessians[g] = ustage->get_constraint_hessians().at(constraint);
                    }
                }
            }
            auto ux = cs::MX::veccat({ret.u, ret.x});
            ret.Gg_eq = generate_jacobian(ux, equality_constraints, jacobians);
            ret.Gg_ineq = generate_jacobian(ux, inequality_constraints, jacobians);
            ret.Gg_dyn = generate_jacobian(ux, x_next, jacobians);
            if (ret.xp1.size1() > 0)
                ret.Gg_dyn.second -= ret.xp1;
            ret.L = 0;
            for (auto &term : ustage->get_objective_terms())
            {
                ret.L += term;
            }
            // add an empty hessian for L
            hessians[ret.L] = Hessian();
            cs::MX lam_dum;
            uo_map_mx<Hessian> hess_dum;
            hess_dum[ret.L] = Hessian();
            ret.hess_obj = generate_hessian(ux, lam_dum, ret.L, {}, hessians);
            ret.hess_dyn = generate_hessian(ux, ret.lam_dyn, cs::MX(0.0), x_next, hessians);
            ret.hess_g_eq = generate_hessian(ux, ret.lam_g_eq, cs::MX(0.0), equality_constraints, hessians);
            ret.hess_g_ineq = generate_hessian(ux, ret.lam_g_ineq, cs::MX(0.0), inequality_constraints, hessians);
            ret.lb = cs::DM::veccat(lb_vec);
            ret.ub = cs::DM::veccat(ub_vec);
            return ret;
        }

        std::pair<cs::MX, cs::MX> uStageQuantities::generate_jacobian(const cs::MX &x, const cs::MX &g, const Jacobian &jac)
        {
            cs::MX ret_G;

            if (jac.is_empty())
            {
                // AD mode
                ret_G = cs::MX::jacobian(g, x);
            }
            else
            {
                if (jac.Jx.size1() != g.size1() || jac.Jx.size2() != x.size1())
                {
                    throw std::runtime_error("Jacobian has wrong dimensions");
                }
                auto P = evalf(jacobian(jac.x, x));
                cs::MX x_r;
                cs::MX jac_r;
                if (P.is_eye())
                {
                    x_r = x;
                    jac_r = jac.Jx;
                }
                else
                {
                    x_r = mtimes(transpose(P), jac.x);
                    jac_r = mtimes(jac.Jx, P);
                }
                auto perm_vec = P.sparsity().permutation_vector(true);
                // check if Jac has the right dimensions
                // check if all elements of x are the same as jac.x
                for (int i = 0; i < jac.x.size1(); i++)
                {
                    if (!is_equal(x(i), jac.x(perm_vec.at(i)), 2))
                    {
                        throw std::runtime_error("x must be the same as in jac");
                    }
                }
                ret_G = jac_r;
            }
            return {ret_G, g};
        }

        std::pair<cs::MX, cs::MX> uStageQuantities::generate_hessian(const cs::MX &x, cs::MX &lam, const cs::MX &J, const cs::MX &g, const Hessian &hess)
        {
            cs::MX ret_H;
            cs::MX ret_h;

            if (hess.is_empty())
            {
                // AD mode
                lam = cs::MX::sym("lam", g.size1(), 1);
                if (g.size1() > 0)
                    ret_H = cs::MX::hessian(J + dot(lam, g), x, ret_h);
                else
                    ret_H = cs::MX::hessian(J, x, ret_h);
            }
            else
            {
                // check if hessian has right dimensions
                if (hess.Hxx.size1() != x.size1() || hess.Hxx.size2() != x.size1())
                {
                    throw std::runtime_error("Hessian has wrong dimensions");
                }
                // x.is_symbolic()
                auto P = evalf(jacobian(hess.x, x));
                cs::MX x_r;
                cs::MX Hx_r;
                cs::MX Hxx_r;
                if (P.is_eye())
                {
                    x_r = x;
                    Hx_r = hess.Hx;
                    Hxx_r = hess.Hxx;
                }
                else
                {
                    x_r = mtimes(transpose(P), hess.x);
                    Hx_r = mtimes(transpose(P), hess.Hx);
                    Hxx_r = mtimes(transpose(P), mtimes(hess.Hxx, P));
                }
                auto perm_vec = P.sparsity().permutation_vector(true);
                // check if all elements of x are the same as hess.x
                for (int i = 0; i < x_r.size1(); i++)
                {
                    if (!is_equal(x(i), hess.x(perm_vec.at(i)), 2))
                    {
                        throw std::runtime_error("x must be the same as in hess");
                    }
                }
                // check if lam has the right dimensions
                if (hess.lam.size() != g.size())
                {
                    throw std::runtime_error("lam has wrong dimensions");
                }
                ret_H = Hxx_r;
                ret_h = Hx_r;
                lam = hess.lam;
            }
            return {ret_H, ret_h};
        }
        std::pair<cs::MX, cs::MX> uStageQuantities::generate_jacobian(const cs::MX &x, const std::vector<cs::MX> &g, const uo_map_mx<Jacobian> &jac)
        {
            const int no_constraints = g.size();
            std::vector<cs::MX> ret_G_vec(no_constraints);
            std::vector<cs::MX> ret_g_vec(no_constraints);
            for (int i = 0; i < no_constraints; i++)
            {
                if (g[i].is_empty())
                    continue;
                auto ret = generate_jacobian(x, g[i], jac.at(g[i]));
                ret_G_vec[i] = ret.first;
                ret_g_vec[i] = ret.second;
            }
            return {cs::MX::vertcat(ret_G_vec), cs::MX::vertcat(ret_g_vec)};
        }
        std::pair<cs::MX, cs::MX> uStageQuantities::generate_hessian(const cs::MX &x, cs::MX &lam, const cs::MX &J, const std::vector<cs::MX> &g, const uo_map_mx<Hessian> &hess)
        {
            const int no_constraints = g.size();
            cs::MX ret_H = cs::MX::zeros(x.size1(), x.size1());
            cs::MX ret_h = cs::MX::zeros(x.size1(), 1);
            std::vector<cs::MX> ret_lam;
            // J
            if (!J.is_zero())
            {
                cs::MX lami;
                auto res = generate_hessian(x, lami, J, {}, hess.at(J));
                ret_H += res.first;
                ret_h += res.second;
            }
            // constraints
            for (int i = 0; i < no_constraints; i++)
            {
                if (g[i].is_empty())
                    continue;
                cs::MX lami;
                auto res = generate_hessian(x, lami, 0.0, g[i], hess.at(g[i]));
                ret_H += res.first;
                ret_h += res.second;
                ret_lam.push_back(lami);
            }
            lam = cs::MX::vertcat(ret_lam);
            return {ret_H, ret_h};
        }
    }
}