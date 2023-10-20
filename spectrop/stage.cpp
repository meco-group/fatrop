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
        }
        void StageInternal::register_control(const cs::MX &control)
        {
            controls_.push_back(control);
            controls_set_.insert(control);
            control_syms_[control] = cs::MX::sym(control.name(), control.size1() * control.size2(), K_);
        }
        void StageInternal::register_control_parameter(const cs::MX &control_parameter)
        {
            control_parameters_.push_back(control_parameter);
            control_parameters_set_.insert(control_parameter);
            control_parameter_syms_[control_parameter] = cs::MX::sym(control_parameter.name(), control_parameter.size1() * control_parameter.size2(), K_);
        }
        bool StageInternal::has_variable(const cs::MX &var)
        {
            bool has_state = states_set_.find(var) != states_set_.end();
            bool has_control = controls_set_.find(var) != controls_set_.end();
            bool has_control_parameter = control_parameters_set_.find(var) != control_parameters_set_.end();
            return has_state || has_control || has_control_parameter;
        }
        const std::vector<cs::MX> &StageInternal::get_objective_terms()
        {
            return objective_terms_;
        };
        const std::vector<cs::MX> &StageInternal::get_constraints()
        {
            return constraints_;
        };
        const uo_map_mx<cs::MX> &StageInternal::get_next_states()
        {
            return next_states_;
        };
        const std::vector<cs::MX> &StageInternal::get_states()
        {
            return states_;
        };
        const uo_map_mx<cs::MX> &StageInternal::get_state_syms()
        {
            return state_syms_;
        };
        const std::vector<cs::MX> &StageInternal::get_controls()
        {
            return controls_;
        };
        const uo_map_mx<cs::MX> &StageInternal::get_control_syms()
        {
            return control_syms_;
        };
        const std::vector<cs::MX> &StageInternal::get_control_parameters()
        {
            return control_parameters_;
        };
        const uo_map_mx<cs::MX> &StageInternal::get_control_parameter_syms()
        {
            return control_parameter_syms_;
        };
        const std::shared_ptr<OcpInternal> &StageInternal::get_ocp()
        {
            return ocp_;
        };
        const std::shared_ptr<StageInternal> &StageInternal::get_next_stage()
        {
            return next_stage_;
        };
    }
}