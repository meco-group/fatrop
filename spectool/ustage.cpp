#include "ustage.hpp"
#include "ocp.hpp"
namespace fatrop
{
    namespace spectool
    {
        void uStage::subject_to(const cs::MX &constraint)
        {
            get()->constraints_.push_back(constraint);
            get()->add_variables(constraint);
        }
        void uStage::add_objective(const cs::MX &objective)
        {
            get()->objective_terms_.push_back(objective);
            get()->add_variables(objective);
        }
        void uStage::set_next(const cs::MX &state, const cs::MX &next_state)
        {
            get()->next_states_[state] = next_state;
            get()->add_variables(next_state);
            if (get()-> next_ustage_.lock())
                get()->next_ustage_.lock()->add_variables(state);
        }
        void uStage::set_next(const uo_map_mx<cs::MX> &next_states)
        {
            for (auto &[state, next_state] : next_states)
                set_next(state, next_state);
        };
        const uo_map_mx<cs::MX> &uStage::dynamics()
        {
            return get()->next_states_;
        }
        void uStageInternal::add_variables(const cs::MX &expr)
        {
            auto syms = cs::MX::symvar(expr);
            for (auto &sym : syms)
            {
                if(!auto_mode)
                {
                    if(!has_variable(sym))
                        throw std::runtime_error("Unregcognized variable in nonauto mode");
                    else
                        continue;
                }
                if (ocp_.lock()->is_global_parameter(sym))
                    continue;
                if (has_variable(sym))
                    continue;
                if (ocp_.lock()->is_state(sym))
                    register_state(sym);
                else if (ocp_.lock()->is_control(sym))
                    register_control(sym);
                else if (ocp_.lock()->is_control_parameter(sym))
                    register_control_parameter(sym);
                else if (ocp_.lock()->is_hybrid(sym))
                    register_hybrid(sym);
                else
                    throw std::runtime_error("MX sym " + sym.name() + " not recognized, is it declared outside of the Ocp?");
            }
        }
        void uStageInternal::register_state(const cs::MX &state)
        {
            states_.push_back(state);
            states_set_.insert(state);
            state_syms_[state] = std::vector<cs::MX>(K_);
            for (int k = 0; k < K_; k++)
                state_syms_[state][k] = cs::MX::sym(state.name() + std::to_string(k), state.size1() * state.size2());
            // state_initial_[state] = cs::DM::zeros(state.size1() * state.size2(), K_);
        }
        void uStageInternal::register_control(const cs::MX &control)
        {
            controls_.push_back(control);
            controls_set_.insert(control);
            control_syms_[control] = std::vector<cs::MX>(K_);
            for (int k = 0; k < K_; k++)
                control_syms_[control][k] = cs::MX::sym(control.name() + std::to_string(k), control.size1() * control.size2());
            // control_initial_[control] = cs::DM::zeros(control.size1() * control.size2(), K_);
        }
        void uStageInternal::register_hybrid(const cs::MX &hybrid)
        {
            hybrids_.push_back(hybrid);
            hybrids_set_.insert(hybrid);
            hybrid_syms_[hybrid] = std::vector<cs::MX>(K_);
            for (int k = 0; k < K_; k++)
                hybrid_syms_[hybrid][k] = cs::MX::sym(hybrid.name() + std::to_string(k), hybrid.size1() * hybrid.size2());
            // automatic_initial_[automatic] = cs::DM::zeros(automatic.size1() * automatic.size2(), K_);
        }
        void uStageInternal::register_control_parameter(const cs::MX &control_parameter)
        {
            control_parameters_.push_back(control_parameter);
            control_parameters_set_.insert(control_parameter);
            control_parameter_syms_[control_parameter] = std::vector<cs::MX>(K_);
            for (int k = 0; k < K_; k++)
                control_parameter_syms_[control_parameter][k] = cs::MX::sym(control_parameter.name() + std::to_string(k), control_parameter.size1() * control_parameter.size2());
            // control_parameter_vals_[control_parameter] = cs::DM::zeros(control_parameter.size1() * control_parameter.size2(), K_);
        }
        void uStageInternal::register_global_parameter(const cs::MX &global_parameter)
        {
            global_parameters_.push_back(global_parameter);
            global_parameters_set_.insert(global_parameter);
        }
        bool uStageInternal::has_variable(const cs::MX &var) const
        {
            bool has_state = states_set_.find(var) != states_set_.end();
            bool has_control = controls_set_.find(var) != controls_set_.end();
            bool has_control_parameter = control_parameters_set_.find(var) != control_parameters_set_.end();
            bool has_hybrid = hybrids_set_.find(var) != hybrids_set_.end();
            // in auto mode global variables are added at the end
            bool has_global = (!auto_mode) && global_parameters_set_.find(var) != global_parameters_set_.end();
            return has_state || has_control || has_control_parameter || has_hybrid || has_global;
        }
        const std::vector<cs::MX> &uStageInternal::get_objective_terms() const
        {
            return objective_terms_;
        };
        const std::vector<cs::MX> &uStageInternal::get_constraints() const
        {
            return constraints_;
        };
        const uo_map_mx<cs::MX> &uStageInternal::get_next_states() const
        {
            return next_states_;
        };
        void uStageInternal::get_hybrids(std::vector<cs::MX> &auto_x, std::vector<cs::MX> &auto_u) const
        {
            for (const auto &hybrid : get_hybrids())
            {
                if (get_prev_ustage())
                {
                    auto &next_states = get_prev_ustage()->get_next_states();
                    (next_states.find(hybrid) != next_states.end() ? auto_x : auto_u).push_back(hybrid);
                }
                else
                {
                    auto_u.push_back(hybrid);
                }
            }
        };
        const std::vector<cs::MX> uStageInternal::get_hybrids_states() const
        {
            std::vector<cs::MX> auto_x;
            std::vector<cs::MX> auto_u;
            get_hybrids(auto_x, auto_u);
            return auto_x;
        };
        const std::vector<cs::MX> uStageInternal::get_hybrids_controls() const
        {
            std::vector<cs::MX> auto_x;
            std::vector<cs::MX> auto_u;
            get_hybrids(auto_x, auto_u);
            return auto_u;
        };
        const std::vector<cs::MX> uStageInternal::get_states(bool include_hybrids) const
        {
            auto ret = states_;
            if (include_hybrids)
            {
                auto auto_x = get_hybrids_states();
                ret.insert(ret.end(), auto_x.begin(), auto_x.end());
            }
            return ret;
        };
        const uo_map_mx<std::vector<cs::MX>> &uStageInternal::get_state_syms() const
        {
            return state_syms_;
        };
        const std::vector<cs::MX> uStageInternal::get_controls(bool include_hybrids) const
        {
            auto ret = controls_;
            if (include_hybrids)
            {
                auto auto_u = get_hybrids_controls();
                ret.insert(ret.end(), auto_u.begin(), auto_u.end());
            }
            return ret;
        };
        const uo_map_mx<std::vector<cs::MX>> &uStageInternal::get_control_syms() const
        {
            return control_syms_;
        };
        const std::vector<cs::MX> &uStageInternal::get_hybrids() const
        {
            return hybrids_;
        };
        const uo_map_mx<std::vector<cs::MX>> &uStageInternal::get_hybrid_syms() const
        {
            return hybrid_syms_;
        };
        const std::vector<cs::MX> &uStageInternal::get_control_parameters() const
        {
            return control_parameters_;
        };
        const uo_map_mx<std::vector<cs::MX>> &uStageInternal::get_control_parameter_syms() const
        {
            return control_parameter_syms_;
        };
        const std::shared_ptr<OcpInternal> uStageInternal::get_ocp() const
        {
            return ocp_.lock();
        };
        const std::shared_ptr<uStageInternal> uStageInternal::get_next_ustage() const
        {
            return next_ustage_.lock();
        };
        const std::shared_ptr<uStageInternal> uStageInternal::get_prev_ustage() const
        {
            return prev_ustage_.lock();
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
        bool uStage::has_variable(const cs::MX &var) const
        {
            return get()->has_variable(var);
        }
        int uStage::K() const
        {
            return get()->K_;
        }
        cs::MX uStage::eval_at_control(const cs::MX &expr, const int k) const
        {
            std::vector<cs::MX> from;
            std::vector<cs::MX> to;
            for (auto &state : get()->get_states(false))
            {
                from.push_back(state);
                to.push_back(get()->get_state_syms().at(state)[k]);
            }
            for (auto &control : get()->get_controls(false))
            {
                from.push_back(control);
                to.push_back(get()->get_control_syms().at(control)[k]);
            }
            for (auto &hybrid : get()->get_hybrids())
            {
                from.push_back(hybrid);
                to.push_back(get()->get_hybrid_syms().at(hybrid)[k]);
            }
            for (auto &control_parameter : get()->get_control_parameters())
            {
                from.push_back(control_parameter);
                to.push_back(get()->get_control_parameter_syms().at(control_parameter)[k]);
            }
            return cs::MX::substitute(std::vector<cs::MX>{expr}, from, to)[0];
        }
        cs::MX uStage::sample(const cs::MX &expr) const
        {
            // check if expr is a column vector
            if (expr.size2() > 1)
                throw std::runtime_error("sample: expr must be a column vector");
            // auto ret = cs::MX::zeros(expr.size1(), K());
            std::vector<cs::MX> samples_vec;
            for (int k = 0; k < K(); k++)
            {
                samples_vec.push_back(eval_at_control(expr, k));
            }
            return cs::MX::horzcat(samples_vec);
        }
        const std::vector<cs::MX> &uStage::get_objective_terms() const
        {
            return get()->get_objective_terms();
        }
        std::shared_ptr<uStageInternal> uStage::get_internal() const
        {
            return static_cast<std::shared_ptr<uStageInternal>>(*this);
        }
        const std::vector<cs::MX> uStage::get_states(bool incl_auto) const
        {
            return get()->get_states(incl_auto);
        }
        const std::vector<cs::MX> uStage::get_controls(bool incl_auto) const
        {
            return get()->get_controls(incl_auto);
        }
        const std::vector<cs::MX> &uStage::get_hybrids() const
        {
            return get()->get_hybrids();
        }
        const std::vector<cs::MX> &uStage::get_control_parameters() const
        {
            return get()->get_control_parameters();
        }
    }
}