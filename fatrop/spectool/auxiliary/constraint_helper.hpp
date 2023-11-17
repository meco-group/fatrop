
#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include <cassert>
#include "casadi_utilities.hpp"
namespace fatrop
{
    namespace spectool
    {
        namespace cs = casadi;
        class ConstraintHelper
        {
        public:
            static void process(const casadi::MX &constr,
                                casadi::DM &lb,
                                casadi::DM &ub,
                                casadi::MX &g_ineq,
                                casadi::MX &g)
            {
                // reset lb, ub, g_ineq, g
                std::vector<cs::DM> lb_vec;
                std::vector<cs::DM> ub_vec;
                std::vector<cs::MX> g_ineq_vec;
                std::vector<cs::MX> g_vec;
                assert(constr.size2() == 1);
                auto opti = casadi::Opti();
                auto syms = casadi::MX::symvar(constr);
                std::vector<casadi::MX> syms_opti;
                for (auto &sym : syms)
                {
                    syms_opti.push_back(opti.variable(sym.size1(), sym.size2()));
                }
                auto canon = opti.advanced().canon_expr(casadi::MX::substitute({constr}, syms, syms_opti)[0]);
                casadi::DM lb_temp = casadi::MX::evalf(canon.lb);
                casadi::DM ub_temp = casadi::MX::evalf(canon.ub);
                auto g_temp = casadi::MX::substitute({canon.canon}, syms_opti, syms)[0];
                for (int i = 0; i < lb_temp.size1(); i++)
                {
                    if ((double)lb_temp(i) == (double)ub_temp(i))
                    {
                        g_vec.push_back(g_temp(i) - lb_temp(i));
                    }
                    else
                    {
                        g_ineq_vec.push_back(g_temp(i));
                        lb_vec.push_back(lb_temp(i));
                        ub_vec.push_back(ub_temp(i));
                    }
                }
                lb = cs::DM::veccat(lb_vec);
                ub = cs::DM::veccat(ub_vec);
                g_ineq = cs::MX::veccat(g_ineq_vec);
                g = cs::MX::veccat(g_vec);
            }
            static void process(const std::vector<casadi::MX> &constr,
                                casadi::DM &lb,
                                casadi::DM &ub,
                                casadi::MX &g_ineq,
                                casadi::MX &g)
            {
                std::vector<cs::DM> lb_vec;
                std::vector<cs::DM> ub_vec;
                std::vector<cs::MX> g_ineq_vec;
                std::vector<cs::MX> g_vec;
                for (auto &constr_i : constr)
                {
                    auto lb_i = casadi::DM();
                    auto ub_i = casadi::DM();
                    auto g_ineq_i = casadi::MX();
                    auto g_i = casadi::MX();
                    process(constr_i, lb_i, ub_i, g_ineq_i, g_i);
                    lb_vec.push_back(lb_i);
                    ub_vec.push_back(ub_i);
                    g_ineq_vec.push_back(g_ineq_i);
                    g_vec.push_back(g_i);
                }
                lb = cs::DM::veccat(lb_vec);
                ub = cs::DM::veccat(ub_vec);
                g_ineq = cs::MX::veccat(g_ineq_vec);
                g = cs::MX::veccat(g_vec);
            }
        };
    }
}