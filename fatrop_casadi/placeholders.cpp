#include "placeholders.hpp"
#include "method-fatrop.hpp"
///////////// IMPLEMENTATION
using namespace fatrop::fatrop_casadi;

std::vector<PlaceHolderType> Placeholders::get_all_types(const casadi::MX &expr)
{
    std::vector<PlaceHolderType> ret;
    for (const auto &p : get_all_placeholders(expr))
    {
        ret.push_back(p.first.type);
    }
    return ret;
}

casadi::MX Placeholders::operator()(const casadi::MX &expr, MXPlaceholder::evaluation_mode mode)
{
    casadi::MX ret = expr;
    while (!get_all_placeholders(ret).empty()) // todo re-use this result
    {
        std::vector<casadi::MX> from;
        std::vector<casadi::MX> to;
        auto placeholders = get_all_placeholders(ret);
        while (!placeholders.empty())
        {
            auto p = placeholders.back();
            placeholders.pop_back();
            from.push_back(p.second);
            to.push_back(p.first.stage->fill_placeholder(p.first.type, p.first, mode));
            // std::cout << "replacing " << p.second << " with " << to.back() << std::endl;
        }
        ret = casadi::MX::substitute({ret}, from, to)[0];
    }
    return ret;
}

std::vector<std::pair<MXPlaceholder, casadi::MX>> Placeholders::get_all_placeholders(const casadi::MX &expr)
{
    auto ret = std::vector<std::pair<MXPlaceholder, casadi::MX>>();
    auto syms = casadi::MX::symvar(expr);
    for (auto &sym : syms)
    {
        if (find(sym) != end())
        {
            ret.push_back(std::make_pair(operator[](sym), sym));
        }
    }
    // std::cout << "number of placeholders: " << ret.size() << std::endl;
    return ret;
}

bool Placeholders::has_placeholders(const casadi::MX &expr, const std::vector<PlaceHolderType> &types)
{
    auto placeholders = get_all_placeholders(expr);
    for (auto &p : placeholders)
    {
        if (std::find(types.begin(), types.end(), p.first.type) != types.end())
            return true;
    }
    return false;
}