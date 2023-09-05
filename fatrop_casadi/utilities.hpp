#pragma once

#include <casadi/casadi.hpp>
namespace fatrop
{
    namespace fatrop_casadi
    {
        struct comp_mx
        {
            bool operator()(const casadi::MX &a, const casadi::MX &b) const
            {
                return a.get() < b.get();
            }
        };
        class ConstraintHelper
        {
        public:
            ConstraintHelper(const casadi::MX &constr)
            {
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
                        g = casadi::MX::veccat({g, g_temp(i)});
                    }
                    else
                    {
                        g_ineq = casadi::MX::veccat({g_ineq, g_temp(i), g_temp(i)});
                        lb = casadi::DM::veccat({lb, lb_temp(i)});
                        lb = casadi::DM::veccat({ub, ub_temp(i)});
                    }
                }
            }
            casadi::DM lb;
            casadi::DM ub;
            casadi::MX g_ineq;
            casadi::MX g;
        };
        class get_terms_helper
        {
        public:
            static std::vector<casadi::MX> get_terms(const casadi::MX &expr)
            {
                std::vector<casadi::MX> ret;
                get_terms_internal(expr, ret);
                return ret;
            }

        private:
            static void get_terms_internal(const casadi::MX &expr, std::vector<casadi::MX> &ret)
            {
                const int n_dep = expr.n_dep();
                if (expr.op() == casadi::OP_ADD)
                {
                    for (int i = 0; i < n_dep; i++)
                    {
                        get_terms_internal(expr.dep(i), ret);
                    }
                }
                else
                {
                    ret.push_back(expr);
                }
            }
        };
    }
    class DM_to_vec_helper
    {
    public:
        static std::vector<double> DM_to_vec(const casadi::DM &dm)
        {
            // check if dm is a vector
            if (dm.size2() != 1)
            {
                throw std::runtime_error("DM_to_vec_helper::DM_to_vec: dm is not a vector");
            }
            std::vector<double> ret(dm.size1());
            for (int i = 0; i < dm.size1(); i++)
            {
                ret.push_back((double)dm(i));
            }
            return ret;
        }
    };
}