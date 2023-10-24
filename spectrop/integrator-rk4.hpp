#pragma once
#include <casadi/casadi.hpp>
#include <vector>
#include <map>
#include "casadi_utilities.hpp"
namespace fatrop
{
    namespace spectrop
    {
        namespace cs = casadi;
        class IntegratorRk4
        {
        public:
            IntegratorRk4(const uo_map_mx<cs::MX>& xdx, const casadi::MX &dt)
            {
                // put the x_dx keys into a vector
                std::vector<cs::MX> x_vec;
                std::vector<cs::MX> dx_vec;
                for (auto &xdx : xdx)
                {
                    x_vec.push_back(xdx.first);
                    dx_vec.push_back(xdx.second);
                }
                x = cs::MX::veccat(x_vec);
                dx = cs::MX::veccat(dx_vec);
                // use RK4 scheme to get x_next
                cs::MX k1 = dt * dx;
                cs::MX k2 = dt * cs::MX::substitute(dx, x, x + 0.5 * k1);
                cs::MX k3 = dt * cs::MX::substitute(dx, x, x + 0.5 * k2);
                cs::MX k4 = dt * cs::MX::substitute(dx, x, x + k3);
                x_next = x + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            }
            cs::MX operator()(const cs::MX &expr)
            {
                return cs::MX::substitute(expr, x, x_next);
            }
            cs::MX dt;
            cs::MX x;
            cs::MX dx;
            cs::MX x_next;
        };
    }
}