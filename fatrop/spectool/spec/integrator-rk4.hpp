#pragma once
#include <casadi/casadi.hpp>
#include <vector>
#include <map>
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
namespace fatrop
{
    namespace spectool
    {
        namespace cs = casadi;
        class IntegratorRk4
        {
        public:
            IntegratorRk4(const std::vector<std::pair<cs::MX, cs::MX>> & xdx_in, const casadi::MX &dt_in)
            {
                dt = dt_in;
                // put the x_dx keys into a vector
                std::vector<cs::MX> x_vec;
                std::vector<cs::MX> dx_vec;
                for (auto &xdx : xdx_in)
                {
                    x_vec.push_back(xdx.first);
                    dx_vec.push_back(xdx.second);
                }
                x = cs::MX::veccat(x_vec);
                dx = cs::MX::veccat(dx_vec);
                // use RK4 scheme to get x_next
                auto f = cs::Function("f", {x}, {dx}, cs::Dict{{"allow_free", true}});
                auto  k1 = dt * f(x)[0];
                auto  k2 = dt * f(x + 0.5 * k1)[0];
                auto  k3 = dt * f(x + 0.5 * k2)[0];
                auto  k4 = dt * f(x + k3)[0];
                // cs::MX k1 = dt * dx;
                // cs::MX k2 = dt * cs::MX::substitute(dx, x, x + 0.5 * k1);
                // cs::MX k3 = dt * cs::MX::substitute(dx, x, x + 0.5 * k2);
                // cs::MX k4 = dt * cs::MX::substitute(dx, x, x + k3);
                x_next = x + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            }
            cs::MX operator()(const cs::MX &expr)
            {
                return operator()(std::vector<cs::MX>{expr})[0];
            }
            std::vector<cs::MX> operator()(const std::vector<cs::MX> &expr)
            {
                if(cs::Function("help", {x}, expr, cs::Dict({{"allow_free", true}})).has_free())
                {
                    throw std::runtime_error("IntegratorRk4: expr has free variables");
                }
                return cs::MX::substitute(expr, {x}, {x_next});
            }

            cs::MX dt;
            cs::MX x;
            cs::MX dx;
            cs::MX x_next;
        };
    }
}