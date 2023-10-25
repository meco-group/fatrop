#include "ocp.hpp"
#include "fatrop_solver.hpp"
namespace fatrop
{
    namespace spectrop
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
        Stage Ocp::new_stage(const int K)
        {
            auto ret = Stage(K, static_cast<std::shared_ptr<OcpInternal>>(*this));
            // add the states to the new stage
            if (!stages_.empty())
            {
                for (auto &[state, sym] : stages_.back().get()->next_states_)
                    ret->add_variables(state);
                // update next stage
                stages_.back().get()->next_stage_ = ret;
            }
            stages_.push_back(ret);
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
            for (const auto &stage : get_stages())
            {
                if (std::all_of(vars.begin(), vars.end(), [&](const cs::MX &var)
                                { return stage.has_variable(var); }))
                    ret = cs::MX::horzcat({ret, stage.sample(expr)});
            }
            return ret;
        }
        const std::vector<Stage> &Ocp::get_stages() const
        {
            return stages_;
        }
        cs::Function Ocp::to_function(const std::vector<cs::MX> &in, const std::vector<cs::MX> &out)
        {
            auto solver = SolverFatrop();
            solver.transcribe(*this);
            std::vector<cs::MX> gist_solver_in;
            std::vector<cs::MX> gist_solver_out;
            auto fatrop_func = solver.to_function(*this, gist_solver_in, gist_solver_out);
            cs::MX vars = gist_solver_in[0];
            auto helper0 = cs::Function("helper0", {vars}, in, cs::Dict{{"allow_free",true}});
            return fatrop_func;
            // if (helper0.has_free())
            // {
            //     auto free_inits = helper0.free_mx();
            //     for (const auto &free_var : free_inits)
            //     {
            //     }
            // }
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