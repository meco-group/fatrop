#pragma once
#include <casadi/casadi.hpp>
#include <string>

namespace fatropy
{
    casadi::Function get_func_from_py(const std::string &file_name, const std::string &func_name);
}