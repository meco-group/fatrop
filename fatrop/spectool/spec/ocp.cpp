#include "ocp.hpp"
#include "fatrop/spectool/solver_interfaces/fatrop/fatrop_solver.hpp"
#include "fatrop/spectool/solver_interfaces/casadi_opti/opti_solver.hpp"
#include <numeric>
namespace fatrop
{
    namespace spectool
    {
        cs::MX Ocp::state(const int m, const int n)
        {
            auto x = cs::MX::sym(std::string("x") + std::to_string(get()->states_.size()), m, n);
            get()->states_.insert(x);
            get()->add_to_ordering(x);
            return x;
        }
        cs::MX Ocp::control(const int m, const int n)
        {
            auto u = cs::MX::sym(std::string("u") + std::to_string(get()->controls_.size()), m, n);
            get()->controls_.insert(u);
            get()->add_to_ordering(u);
            return u;
        }
        cs::MX Ocp::hybrid(const int m, const int n)
        {
            auto hybrid = cs::MX::sym(std::string("hybrid") + std::to_string(get()->hybrids_.size()), m, n);
            get()->hybrids_.insert(hybrid);
            get()->add_to_ordering(hybrid);
            return hybrid;
        }
        cs::MX Ocp::parameter(const int m, const int n, const std::string &grid)
        {
            if (grid == "global")
            {
                auto p = cs::MX::sym(std::string("p") + std::to_string(get()->global_parameters_.size()), m, n);
                get()->global_parameters_.insert(p);
                get()->global_parammeter_syms_.push_back(p);
                get()->add_to_ordering(p);
                return p;
            }
            if (grid == "control")
            {
                auto p = cs::MX::sym(std::string("p") + std::to_string(get()->control_parameters_.size()), m, n);
                get()->control_parameters_.insert(p);
                get()->add_to_ordering(p);
                return p;
            }
            throw std::runtime_error("grid must be either 'global' or 'control'");
            return cs::MX();
        }
        uStage Ocp::new_ustage(const int K)
        {
            auto ret = uStage(K, static_cast<std::shared_ptr<OcpInternal>>(*this));
            // add the states to the new stage
            ustages_.push_back(ret);
            return ret;
        }

        void Ocp::add_ustage(const uStage &ustage)
        {
            // check if ustage is already registered
            if (std::find(ustages_.begin(), ustages_.end(), ustage) != ustages_.end())
                throw std::runtime_error("ustage is already registered, consider using uStage::clone()");
            ustages_.push_back(ustage);
        }

        bool OcpInternal::is_state(const cs::MX &var)
        {
            return states_.find(var) != states_.end();
        }
        bool OcpInternal::is_control(const cs::MX &var)
        {
            return controls_.find(var) != controls_.end();
        }
        bool OcpInternal::is_hybrid(const cs::MX &var)
        {
            return hybrids_.find(var) != hybrids_.end();
        }
        bool OcpInternal::is_global_parameter(const cs::MX &var)
        {
            return global_parameters_.find(var) != global_parameters_.end();
        }
        bool OcpInternal::is_control_parameter(const cs::MX &var)
        {
            return control_parameters_.find(var) != control_parameters_.end();
        }
        const uo_set_mx &OcpInternal::get_states()
        {
            return states_;
        }
        const uo_set_mx &OcpInternal::get_controls()
        {
            return controls_;
        }
        const uo_set_mx &OcpInternal::get_global_parameters()
        {
            return global_parameters_;
        }
        const uo_set_mx &OcpInternal::get_control_parameters()
        {
            return control_parameters_;
        }
        void OcpInternal::add_to_ordering(const cs::MX &var)
        {
            ordering_[var] = ordering_.size();
        }
        std::vector<cs::MX> OcpInternal::order_vars(const std::vector<cs::MX> &vars)
        {
            auto ret = vars;
            std::sort(ret.begin(), ret.end(), [this](const cs::MX &a, const cs::MX &b)
                      { return ordering_[a] < ordering_[b]; });
            return ret;
        }
        std::pair<std::vector<int>, cs::MX> Ocp::sample(const cs::MX &expr) const
        {
            if (expr.size2() != 1)
                return {{}, cs::MX()}; // return empty matrix if input is not a column vector
            auto ret = std::vector<cs::MX>();
            auto vars = cs::MX::symvar(expr);
            auto reti = std::vector<int>();
            int k = 0;
            for (size_t i = 0; i < get_ustages().size(); ++i)
            {
                const auto &ustage = get_ustages()[i];
                {
                    auto sample_ = ustage.sample(expr);
                    if (sample_.second.size2() == 0)
                        continue;
                    // add samples to ret
                    ret.push_back(sample_.second);
                    // add indices to reti
                    for (const auto &idx : sample_.first)
                        reti.push_back(k + idx);
                }
                k += ustage.K();
            }
            return {reti, cs::MX::horzcat(ret)};
        }
        const std::vector<uStage> &Ocp::get_ustages() const
        {
            return ustages_;
        }
        cs::Function Ocp::to_function(const std::string &name, const std::vector<cs::MX> &in, const std::vector<cs::MX> &out) const
        {
            if (get()->solver_name == "")
                throw std::runtime_error("solver not set, use ocp.solver(\"fatrop\") or ocp.solver(\"ipopt\")");
            if (get()->solver_name == "fatrop")
                get()->solver_ptr = std::make_shared<SolverFatrop>();
            else if (get()->solver_name == "ipopt")
                get()->solver_ptr = std::make_shared<SolverOpti>();
            else
                throw std::runtime_error("solver not supported");
            get()->solver_ptr->transcribe(*this, function_opts_);
            std::vector<cs::MX> gist_solver_in;
            std::vector<cs::MX> gist_solver_out;
            auto fatrop_func = get()->solver_ptr->to_function(name, *this, gist_solver_in, gist_solver_out, solver_opts_);
            cs::MX vars = gist_solver_in[0];
            cs::MX initial_guess = gist_solver_in[0];
            auto helper0 = cs::Function("helper0", in, {vars}, cs::Dict{{"allow_free", true}});
            auto helper1 = cs::Function("helper1", in, {gist_solver_in[1], gist_solver_in[2]}, cs::Dict{{"allow_free", true}});
            if (helper1.has_free())
                throw std::runtime_error("to_function: problem still has undefined parameters");
            if (helper0.has_free())
            {
                auto free_inits = helper0.free_mx();
                std::vector<cs::MX> init_evals = eval_at_initial(free_inits);
                initial_guess = cs::MX::substitute({initial_guess}, free_inits, init_evals)[0];
            }
            auto result = fatrop_func({initial_guess, gist_solver_in[1], gist_solver_in[2]});
            return cs::Function(name, in, cs::MX::substitute(out, gist_solver_out, result));
        }
        void Ocp::set_initial(const cs::MX &var, const cs::MX &value)
        {
            get()->initial_values.push_back({var, value});
        }
        std::vector<cs::MX> Ocp::eval_at_initial(const std::vector<cs::MX> &expr) const
        {
            auto all_vars_ocp = all_variables();
            auto all_vars = cs::MX::symvar(veccat(expr));
            std::vector<cs::MX> from;
            std::vector<cs::MX> to;
            for (const auto &[var, value] : get()->initial_values)
            {
                cs::MX varr;
                cs::MX valuee;
                if (var.size1() != value.size1())
                    throw std::runtime_error("initial value has wrong size");

                if (get()->is_state(var) || get()->is_control(var) || get()->is_hybrid(var) || get()->is_control_parameter(var))
                    varr = sample(var).second;
                else if (get()->is_global_parameter(var))
                    varr = var;
                else if (!cs::Function("helper", {all_vars_ocp}, {var}, cs::Dict{{"allow_free", true}}).has_free())
                    varr = var;
                else
                    throw std::runtime_error("unrecognized variable at eval_at_initial");
                if (value.size2() == 1)
                    valuee = cs::MX::repmat(value, 1, varr.size2());
                else
                    valuee = value;
                from.push_back(varr);
                to.push_back(valuee);
            }
            // set the values for which no initialization is provided to zero
            auto helper0 = cs::Function("helper0", from, all_vars, cs::Dict{{"allow_free", true}});
            for (const auto &free : helper0.free_mx())
            {
                from.push_back(free);
                to.push_back(cs::DM::zeros(free.size1(), free.size2()));
            }

            return cs::MX::substitute(expr, from, to);
        }
        Stage Ocp::new_stage(const int K)
        {
            return Stage(*this, K);
        }
        cs::MX Ocp::all_variables() const
        {
            std::vector<cs::MX> gist_in;
            std::vector<cs::MX> gist_out;
            SolverFatrop().gist(*this, gist_in, gist_out);
            return cs::MX::veccat(gist_out);
        }

        void Ocp::subject_to(const cs::MX &expr)
        {
            auto vars = cs::MX::symvar(expr);
            std::vector<uStage> dep_ustage_list;
            // iterate over all uStages and check if they depend on the variables in expr
            for (const auto &ustage : ustages_)
            {
                // auto ustage_vars = ustage.get_variables();
                if (ustage.K() == 1 && ustage.get()->is_evaluable(expr))
                    dep_ustage_list.push_back(ustage);
            }
            // if only one ustage
            if (dep_ustage_list.size() == 1)
            // add to this ustage
            {
                auto to = cs::MX::veccat(dep_ustage_list[0].all_vars());
                auto from = dep_ustage_list[0].at_t0(to);
                dep_ustage_list[0].subject_to(cs::MX::substitute(expr, from, to));
            }
            else
            {
                // if multiple uStages add to multi_ustage_constraints
                get()->multi_ustage_constraints.push_back(expr);
            }
        }
        // void Ocp::set_initial(const cs::MX &var, const cs::MX &initial)
        // {
        //     bool column_mode = initial.size2() == 1;
        //     int offs = 0;
        //     // iterate over stages and set initial
        //     for (auto &stage : stages_)
        //     {
        //         if (stage.has_variable(var))
        //         {
        //             if (column_mode)
        //                 stage.set_initial(var, var);
        //             else
        //             {
        //                 if (offs + stage.K() > initial.size2())
        //                     throw std::runtime_error("initial value has wrong size");
        //                 stage.set_initial(var, initial(cs::Slice(), cs::Slice(offs, offs + stage.K())));
        //             }
        //             offs += stage.K();
        //         }
        //     }
        //     if (!column_mode && offs != initial.size2())
        //         throw std::runtime_error("initial value has wrong size");
        // }
        // void Ocp::set_value(const cs::MX &var, const cs::DM &val)
        // {
        //     bool column_mode = val.size2() == 1;
        //     int offs = 0;
        //     // iterate over stages and set value
        //     for (auto &stage : stages_)
        //     {
        //         if (stage.has_variable(var))
        //         {
        //             if (column_mode)
        //                 stage.set_value(var, val);
        //             else
        //             {
        //                 if (offs + stage.K() > val.size2())
        //                     throw std::runtime_error("initial value has wrong size");
        //                 stage.set_initial(var, val(cs::Slice(), cs::Slice(offs, offs + stage.K())));
        //             }
        //             offs += stage.K();
        //         }
        //     }
        //     if (!column_mode && offs != val.size2())
        //         throw std::runtime_error("value has wrong size");
        // }
    } // namespace spectrop
} // namespace fatrop