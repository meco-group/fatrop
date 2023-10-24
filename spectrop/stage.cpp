#include "stage.hpp"
#include "ocp.hpp"
namespace fatrop
{
    namespace spectrop
    {
        void Stage::subject_to(const cs::MX &constraint)
        {
            get()->constraints_.push_back(constraint);
            get()->add_variables(constraint);
        }
        void Stage::add_objective(const cs::MX &objective)
        {
            get()->objective_terms_.push_back(objective);
            get()->add_variables(objective);
        }
        void Stage::set_next(const cs::MX &state, const cs::MX &next_state)
        {
            get()->next_states_[state] = next_state;
            get()->add_variables(next_state);
            if (get()->next_stage_)
                get()->next_stage_->add_variables(state);
        }
        void Stage::set_next(const uo_map_mx<cs::MX> &next_states)
        {
            for (auto &[state, next_state] : next_states)
                set_next(state, next_state);
        };
        const uo_map_mx<cs::MX> &Stage::dynamics()
        {
            return get()->next_states_;
        }
        void StageInternal::add_variables(const cs::MX &expr)
        {
            auto syms = cs::MX::symvar(expr);
            for (auto &sym : syms)
            {
                if (ocp_->is_global_parameter(sym))
                    continue;
                if (has_variable(sym))
                    continue;
                if (ocp_->is_state(sym))
                    register_state(sym);
                else if (ocp_->is_control(sym))
                    register_control(sym);
                else if (ocp_->is_control_parameter(sym))
                    register_control_parameter(sym);
                else
                    throw std::runtime_error("MX sym " + sym.name() + " not recognized, is it declared outside of the Ocp?");
            }
        }
        void StageInternal::register_state(const cs::MX &state)
        {
            states_.push_back(state);
            states_set_.insert(state);
            state_syms_[state] = cs::MX::sym(state.name(), state.size1() * state.size2(), K_);
            // state_initial_[state] = cs::DM::zeros(state.size1() * state.size2(), K_);
        }
        void StageInternal::register_control(const cs::MX &control)
        {
            controls_.push_back(control);
            controls_set_.insert(control);
            control_syms_[control] = cs::MX::sym(control.name(), control.size1() * control.size2(), K_);
            // control_initial_[control] = cs::DM::zeros(control.size1() * control.size2(), K_);
        }
        void StageInternal::register_control_parameter(const cs::MX &control_parameter)
        {
            control_parameters_.push_back(control_parameter);
            control_parameters_set_.insert(control_parameter);
            control_parameter_syms_[control_parameter] = cs::MX::sym(control_parameter.name(), control_parameter.size1() * control_parameter.size2(), K_);
            // control_parameter_vals_[control_parameter] = cs::DM::zeros(control_parameter.size1() * control_parameter.size2(), K_);
        }
        bool StageInternal::has_variable(const cs::MX &var) const
        {
            bool has_state = states_set_.find(var) != states_set_.end();
            bool has_control = controls_set_.find(var) != controls_set_.end();
            bool has_control_parameter = control_parameters_set_.find(var) != control_parameters_set_.end();
            return has_state || has_control || has_control_parameter;
        }
        const std::vector<cs::MX> &StageInternal::get_objective_terms() const
        {
            return objective_terms_;
        };
        const std::vector<cs::MX> &StageInternal::get_constraints() const
        {
            return constraints_;
        };
        const uo_map_mx<cs::MX> &StageInternal::get_next_states() const
        {
            return next_states_;
        };
        const std::vector<cs::MX> &StageInternal::get_states() const
        {
            return states_;
        };
        const uo_map_mx<cs::MX> &StageInternal::get_state_syms() const
        {
            return state_syms_;
        };
        const std::vector<cs::MX> &StageInternal::get_controls() const
        {
            return controls_;
        };
        const uo_map_mx<cs::MX> &StageInternal::get_control_syms() const
        {
            return control_syms_;
        };
        const std::vector<cs::MX> &StageInternal::get_control_parameters() const
        {
            return control_parameters_;
        };
        const uo_map_mx<cs::MX> &StageInternal::get_control_parameter_syms() const
        {
            return control_parameter_syms_;
        };
        const std::shared_ptr<OcpInternal> &StageInternal::get_ocp() const
        {
            return ocp_;
        };
        const std::shared_ptr<StageInternal> &StageInternal::get_next_stage() const
        {
            return next_stage_;
        };
        // const uo_map_mx<cs::MX> &StageInternal::get_state_initial() const
        // {
        //     return state_initial_;
        // }
        // const uo_map_mx<cs::MX> &StageInternal::get_control_initial() const
        // {
        //     return control_initial_;
        // }
        // const uo_map_mx<cs::DM> &StageInternal::get_control_parameter_vals() const
        // {
        //     return control_parameter_vals_;
        // }
        // void Stage::set_initial(const cs::MX &var, const cs::MX &initial)
        // {
        //     // check if initial is the right size
        //     if (initial.size1() != var.size1() || initial.size2() != get()->K_)
        //         throw std::runtime_error("initial value for state " + var.name() + " has the wrong size");
        //     // check if var is a state or control
        //     if (get()->states_set_.find(var) != get()->states_set_.end())
        //     {
        //         get()->state_initial_[var] = initial;
        //         return;
        //     }
        //     if (get()->controls_set_.find(var) != get()->controls_set_.end())
        //     {
        //         get()->control_initial_[var] = initial;
        //         return;
        //     }
        //     throw std::runtime_error("MX sym " + var.name() + " not recognized as a control or state, is it declared outside of the Ocp?");
        // }
        // void Stage::set_value(const cs::MX &var, const cs::DM &val)
        // {
        //     if (val.size2() == 1 && get()->K_ != 1) // if val is a vector, make it a matrix
        //         return set_value(var, cs::DM::repmat(val, 1, get()->K_));
        //     // check if initial is the right size
        //     if (val.size1() != var.size1() || val.size2() != get()->K_)
        //         throw std::runtime_error("initial value for state " + var.name() + " has the wrong size");
        //     if (get()->control_parameter_vals_.find(var) != get()->control_parameter_vals_.end())
        //     {
        //         get()->control_parameter_vals_[var] = val;
        //         return;
        //     }
        //     throw std::runtime_error("MX sym " + var.name() + " not recognized as a control-grid parameter, is it declared outside of the Ocp?");
        // }
        bool Stage::has_variable(const cs::MX &var) const
        {
            return get()->has_variable(var);
        }
        int Stage::K() const
        {
            return get()->K_;
        }
        cs::MX Stage::eval_at_control(const cs::MX &expr, const int k) const
        {
            std::vector<cs::MX> from;
            std::vector<cs::MX> to;
            for (auto &state : get()->get_states())
            {
                from.push_back(state);
                to.push_back(state(cs::Slice(), k));
            }
            for (auto &control : get()->get_controls())
            {
                from.push_back(control);
                to.push_back(control(cs::Slice(), k));
            }
            for (auto &control_parameter : get()->get_control_parameters())
            {
                from.push_back(control_parameter);
                to.push_back(control_parameter(cs::Slice(), k));
            }
            return cs::MX::substitute(std::vector<cs::MX>{expr}, from, to)[0];
        }
        cs::MX Stage::sample(const cs::MX &expr) const
        {
            // check if expr is a column vector
            if (expr.size2() > 1)
                throw std::runtime_error("sample: expr must be a column vector");
            auto ret = cs::MX::zeros(expr.size1(), K());
            for (int k = 0; k < K(); k++)
            {
                ret(cs::Slice(), k) = eval_at_control(expr, k);
            }
            return ret;
        }
    }
}