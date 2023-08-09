#pragma once
#include <casadi/casadi.hpp>
#include <vector>
#include <map>
#include "auxiliary.hpp"
namespace fatrop_casadi
{
    class IntegratorRk4
    {
    public:
        IntegratorRk4(const std::map<casadi::MX, casadi::MX, comp_mx> &xdx, const casadi::MX &dt)
        {
            // put the x_dx keys into a vector
            std::vector<casadi::MX> x_vec;
            std::vector<casadi::MX> dx_vec;
            for (auto &xdx : xdx)
            {
                x_vec.push_back(xdx.first);
                dx_vec.push_back(xdx.second);
            }
            x = casadi::MX::veccat(x_vec);
            dx = casadi::MX::veccat(dx_vec);
            // use RK4 scheme to get x_next
            casadi::MX k1 = dt * dx;
            casadi::MX k2 = dt * casadi::MX::substitute(dx, x, x + 0.5 * k1);
            casadi::MX k3 = dt * casadi::MX::substitute(dx, x, x + 0.5 * k2);
            casadi::MX k4 = dt * casadi::MX::substitute(dx, x, x + k3);
            x_next = x + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        }
        casadi::MX operator()(const casadi::MX &expr)
        {
            return casadi::MX::substitute(expr, x, x_next);
        }
        casadi::MX dt;
        casadi::MX x;
        casadi::MX dx;
        casadi::MX x_next;
    };
}