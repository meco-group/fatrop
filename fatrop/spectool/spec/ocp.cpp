#include "ocp.hpp"
#include "fatrop/spectool/solver_interfaces/fatrop/fatrop_solver.hpp"
namespace fatrop
{
    namespace spectool
    {
        cs::MX Ocp::state(const int m, const int n)
        {
            auto x = cs::MX::sym(std::string("x") + std::to_string(get()->states_.size()), m, n);
            get()->states_.insert(x);
            return x;
        }
        cs::MX Ocp::control(const int m, const int n)
        {
            auto u = cs::MX::sym(std::string("u") + std::to_string(get()->controls_.size()), m, n);
            get()->controls_.insert(u);
            return u;
        }
        cs::MX Ocp::hybrid(const int m, const int n)
        {
            auto hybrid = cs::MX::sym(std::string("hybrid") + std::to_string(get()->hybrids_.size()), m, n);
            get()->hybrids_.insert(hybrid);
            return hybrid;
        }
        cs::MX Ocp::parameter(const int m, const int n, const std::string &grid)
        {
            if (grid == "global")
            {
                auto p = cs::MX::sym(std::string("p") + std::to_string(get()->global_parameters_.size()), m, n);
                get()->global_parameters_.insert(p);
                get()->global_parammeter_syms_.push_back(p);
                return p;
            }
            if (grid == "control")
            {
                auto p = cs::MX::sym(std::string("p") + std::to_string(get()->control_parameters_.size()), m, n);
                get()->control_parameters_.insert(p);
                return p;
            }
            throw std::runtime_error("grid must be either 'global' or 'control'");
            return cs::MX();
        }
        uStage Ocp::new_ustage(const int K)
        {
            auto ret = uStage(K, static_cast<std::shared_ptr<OcpInternal>>(*this));
            // add the states to the new stage
            if (!ustages_.empty())
            {
                for (auto &[state, sym] : ustages_.back().get()->next_states_)
                    ret->add_variables(state);
                ret->prev_ustage_ = ustages_.back();
                // update next stage
                ustages_.back().get()->next_ustage_ = ret;
            }
            ustages_.push_back(ret);
            return ret;
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
        cs::MX Ocp::sample(const cs::MX &expr) const
        {
            if (expr.size2() != 1)
                return cs::MX(); // return empty matrix if input is not a column vector
            auto ret = cs::MX::zeros(expr.size1(), 0);
            auto vars = cs::MX::symvar(expr);
            for (const auto &ustage : get_ustages())
            {
                if (std::all_of(vars.begin(), vars.end(), [&](const cs::MX &var)
                                { return get()->is_global_parameter(var) || ustage.has_variable(var); }))
                    ret = cs::MX::horzcat({ret, ustage.sample(expr)});
            }
            return ret;
        }
        const std::vector<uStage> &Ocp::get_ustages() const
        {
            return ustages_;
        }
        cs::Function Ocp::to_function(const std::vector<cs::MX> &in, const std::vector<cs::MX> &out, const cs::Dict &opts) const
        {
            auto solver = SolverFatrop();
            solver.transcribe(*this, opts);
            std::vector<cs::MX> gist_solver_in;
            std::vector<cs::MX> gist_solver_out;
            auto fatrop_func = solver.to_function(*this, gist_solver_in, gist_solver_out);
            cs::MX vars = gist_solver_in[0];
            cs::MX initial_guess = eval_at_initial(gist_solver_in[0]);
            auto helper0 = cs::Function("helper0", in, {vars}, cs::Dict{{"allow_free", true}});
            auto helper1 = cs::Function("helper0", in, {gist_solver_in[1], gist_solver_in[2]}, cs::Dict{{"allow_free", true}});
            if (helper1.has_free())
                throw std::runtime_error("to_function: problem still has undefined parameters");
            if (helper0.has_free())
            {
                auto free_inits = helper0.free_mx();
                std::vector<cs::MX> free_zeros;
                for (const auto &free_init : free_inits)
                    free_zeros.push_back(cs::DM::zeros(free_init.size1(), free_init.size2()));
                initial_guess = cs::MX::substitute({initial_guess}, free_inits, free_zeros)[0];
            }
            auto result = fatrop_func({initial_guess, gist_solver_in[1], gist_solver_in[2]});
            return cs::Function("ocp", in, cs::MX::substitute(out, gist_solver_out, result));
        }
        void Ocp::set_initial(const cs::MX &var, const cs::MX &value)
        {
            get()->initial_values.push_back({var, value});
        }
        cs::MX Ocp::eval_at_initial(const cs::MX &expr) const
        {
            std::vector<cs::MX> from;
            std::vector<cs::MX> to;
            for (const auto &[var, value] : get()->initial_values)
            {
                cs::MX varr;
                cs::MX valuee;
                if (var.size1() != value.size1())
                    throw std::runtime_error("initial value has wrong size");

                if (get()->is_state(var) || get()->is_control(var) || get()->is_hybrid(var))
                    varr = sample(var);
                else
                    varr = var;
                if (value.size2() == 1)
                    valuee = cs::MX::repmat(value, 1, varr.size2());
                else
                    valuee = value;
                from.push_back(varr);
                to.push_back(valuee);
            }
            return cs::MX::substitute({expr}, from, to)[0];
        }
        Stage Ocp::new_stage(const int K)
        {
            return Stage(*this, K);
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