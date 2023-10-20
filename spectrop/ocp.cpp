#include "ocp.hpp"
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
    } // namespace spectrop
} // namespace fatrop