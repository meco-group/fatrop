#include "ustage.hpp"
#include "ocp.hpp"
namespace fatrop
{
    namespace spectool
    {
        void uStage::subject_to(const cs::MX &constraint, const Jacobian &jac, const Hessian &hess)
        {
            get()->set_dirty();
            get()->constraints_.push_back(constraint);
            get()->constraint_jacobians_[constraint] = jac;
            get()->constraint_hessians_[constraint] = hess;
            get()->add_variables(constraint);
        }
        void uStage::add_objective(const cs::MX &objective)
        {
            get()->set_dirty();
            get()->objective_terms_.push_back(objective);
            get()->add_variables(objective);
        }
        void uStage::set_next(const cs::MX &state, const cs::MX &next_state, const Jacobian &jac, const Hessian &hess)
        {
            get()->set_dirty();
            get()->next_states_[state] = next_state;
            get()->next_state_jacobians_[state] = jac;
            get()->next_state_hessians_[state] = hess;
            if (K() > 1)
                get()->add_variables(state);
            get()->add_variables(next_state);
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
                if (!auto_mode)
                {
                    if (!has_variable(sym))
                        throw std::runtime_error("Unregcognized variable in nonauto mode");
                    else
                        continue;
                }
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
                else if (ocp_.lock()->is_global_parameter(sym))
                    register_global_parameter(sym);
                else
                    throw std::runtime_error("MX sym " + sym.name() + " not recognized, is it declared outside of the Ocp?");
            }
        }
        void uStageInternal::register_state(const cs::MX &state)
        {
            // check if already registered
            if (states_set_.find(state) == states_set_.end())
            {
                states_.push_back(state);
                if (ocp_.lock())
                    states_ = ocp_.lock()->order_vars(states_);
                states_set_.insert(state);
                set_dirty();
            }
            // check if syms already registered
            if (state_syms_.find(state) == state_syms_.end())
            {
                state_syms_[state] = std::vector<cs::MX>(K_);
                for (int k = 0; k < K_; k++)
                {
                    state_syms_[state][k] = cs::MX::sym(state.name() + std::to_string(k), state.size1() * state.size2());
                    all_eval_syms_.push_back(state_syms_[state][k]);
                }
            }
            // state_initial_[state] = cs::DM::zeros(state.size1() * state.size2(), K_);
        }
        void uStageInternal::register_control(const cs::MX &control)
        {
            // check if already registered
            if (controls_set_.find(control) == controls_set_.end())
            {
                controls_.push_back(control);
                if (ocp_.lock())
                    controls_ = ocp_.lock()->order_vars(controls_);
                controls_set_.insert(control);
                set_dirty();
            }
            // check if syms already registered
            if (control_syms_.find(control) == control_syms_.end())
            {
                control_syms_[control] = std::vector<cs::MX>(K_);
                for (int k = 0; k < K_; k++)
                {
                    control_syms_[control][k] = cs::MX::sym(control.name() + std::to_string(k), control.size1() * control.size2());
                    all_eval_syms_.push_back(control_syms_[control][k]);
                }
            }
            // control_initial_[control] = cs::DM::zeros(control.size1() * control.size2(), K_);
        }
        void uStageInternal::register_hybrid(const cs::MX &hybrid)
        {
            // check if already registered
            if (hybrids_set_.find(hybrid) == hybrids_set_.end())
            {
                hybrids_.push_back(hybrid);
                if (ocp_.lock())
                    hybrids_ = ocp_.lock()->order_vars(hybrids_);
                hybrids_set_.insert(hybrid);
                set_dirty();
            }
            // check if syms already registered
            if (hybrid_syms_.find(hybrid) == hybrid_syms_.end())
            {
                hybrid_syms_[hybrid] = std::vector<cs::MX>(K_);
                for (int k = 0; k < K_; k++)
                {
                    hybrid_syms_[hybrid][k] = cs::MX::sym(hybrid.name() + std::to_string(k), hybrid.size1() * hybrid.size2());
                    all_eval_syms_.push_back(hybrid_syms_[hybrid][k]);
                }
            }
            has_hybrids = true;
            // automatic_initial_[automatic] = cs::DM::zeros(automatic.size1() * automatic.size2(), K_);
        }
        void uStageInternal::register_control_parameter(const cs::MX &control_parameter)
        {
            // check if already registered
            if (control_parameters_set_.find(control_parameter) == control_parameters_set_.end())
            {
                control_parameters_.push_back(control_parameter);
                if (ocp_.lock())
                    control_parameters_ = ocp_.lock()->order_vars(control_parameters_);
                control_parameters_set_.insert(control_parameter);
                set_dirty();
            }
            // check if syms already registered
            if (control_parameter_syms_.find(control_parameter) == control_parameter_syms_.end())
            {
                control_parameter_syms_[control_parameter] = std::vector<cs::MX>(K_);
                for (int k = 0; k < K_; k++)
                {
                    control_parameter_syms_[control_parameter][k] = cs::MX::sym(control_parameter.name() + std::to_string(k), control_parameter.size1() * control_parameter.size2());
                    all_eval_syms_.push_back(control_parameter_syms_[control_parameter][k]);
                }
            }
            // control_parameter_vals_[control_parameter] = cs::DM::zeros(control_parameter.size1() * control_parameter.size2(), K_);
        }
        void uStageInternal::register_global_parameter(const cs::MX &global_parameter)
        {
            // check if already registered
            if (global_parameters_set_.find(global_parameter) == global_parameters_set_.end())
            {
                global_parameters_.push_back(global_parameter);
                if (ocp_.lock())
                    global_parameters_ = ocp_.lock()->order_vars(global_parameters_);
                global_parameters_set_.insert(global_parameter);
                all_eval_syms_.push_back(global_parameter);
            }
        }
        void uStageInternal::register_state(const std::vector<cs::MX> &states)
        {
            std::for_each(states.begin(), states.end(), [&](const auto &state)
                          { register_state(state); });
        }
        void uStageInternal::register_control(const std::vector<cs::MX> &controls)
        {
            std::for_each(controls.begin(), controls.end(), [&](const auto &control)
                          { register_control(control); });
        }
        void uStageInternal::register_hybrid(const std::vector<cs::MX> &hybrids)
        {
            std::for_each(hybrids.begin(), hybrids.end(), [&](const auto &hybrid)
                          { register_hybrid(hybrid); });
        }
        void uStageInternal::register_control_parameter(const std::vector<cs::MX> &control_parameters)
        {
            std::for_each(control_parameters.begin(), control_parameters.end(), [&](const auto &control_parameter)
                          { register_control_parameter(control_parameter); });
        }
        void uStageInternal::register_global_parameter(const std::vector<cs::MX> &global_parameters)
        {
            std::for_each(global_parameters.begin(), global_parameters.end(), [&](const auto &global_parameter)
                          { register_global_parameter(global_parameter); });
        }
        bool uStageInternal::has_variable(const cs::MX &var) const
        {
            bool has_state = states_set_.find(var) != states_set_.end();
            bool has_control = controls_set_.find(var) != controls_set_.end();
            bool has_control_parameter = control_parameters_set_.find(var) != control_parameters_set_.end();
            bool has_hybrid = hybrids_set_.find(var) != hybrids_set_.end();
            bool has_global = global_parameters_set_.find(var) != global_parameters_set_.end();
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
        void uStageInternal::get_hybrids(std::vector<cs::MX> &auto_x, std::vector<cs::MX> &auto_u, const std::shared_ptr<const uStageInternal> &prev) const
        {
            if (has_hybrids && !prev)
                throw std::runtime_error("get_hybrids: prev ustage must be provided if stage has hybrids");

            for (const auto &hybrid : get_hybrids())
            {
                // a hybrid is a state if it is defined by the dynamics function of the previous ustage OR if K>1 and it is defined by the dynamics function of the current ustage
                auto &next_states = prev->get_next_states();
                auto &curr_next_states = get_next_states();
                ((prev && next_states.find(hybrid) != next_states.end()) || (K() > 1 && curr_next_states.find(hybrid) != curr_next_states.end()) ? auto_x : auto_u).push_back(hybrid);
            }
        };
        const std::vector<cs::MX> uStageInternal::get_hybrids_states(const std::shared_ptr<const uStageInternal> &prev) const
        {
            std::vector<cs::MX> auto_x;
            std::vector<cs::MX> auto_u;
            get_hybrids(auto_x, auto_u, prev);
            return auto_x;
        };
        const std::vector<cs::MX> uStageInternal::get_hybrids_controls(const std::shared_ptr<const uStageInternal> &prev) const
        {
            std::vector<cs::MX> auto_x;
            std::vector<cs::MX> auto_u;
            get_hybrids(auto_x, auto_u, prev);
            return auto_u;
        };
        const std::vector<cs::MX> uStageInternal::get_states(bool include_hybrids, const std::shared_ptr<const uStageInternal> &prev) const
        {
            auto ret = states_;
            if (include_hybrids)
            {
                auto auto_x = get_hybrids_states(prev);
                ret.insert(ret.end(), auto_x.begin(), auto_x.end());
            }
            return ret;
        };
        const uo_map_mx<std::vector<cs::MX>> &uStageInternal::get_state_syms() const
        {
            return state_syms_;
        };
        const std::vector<cs::MX> uStageInternal::get_controls(bool include_hybrids, const std::shared_ptr<const uStageInternal> &prev) const
        {
            auto ret = controls_;
            if (include_hybrids)
            {
                auto auto_u = get_hybrids_controls(prev);
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
        const std::vector<cs::MX> &uStageInternal::get_global_parameters() const
        {
            return global_parameters_;
        };
        const uo_map_mx<std::vector<cs::MX>> &uStageInternal::get_control_parameter_syms() const
        {
            return control_parameter_syms_;
        };
        const std::shared_ptr<OcpInternal> uStageInternal::get_ocp() const
        {
            return ocp_.lock();
        };
        const uo_map_mx<Jacobian> &uStageInternal::get_constraint_jacobians() const
        {
            return constraint_jacobians_;
        }
        const uo_map_mx<Hessian> &uStageInternal::get_constraint_hessians() const
        {
            return constraint_hessians_;
        }
        const uo_map_mx<Jacobian> &uStageInternal::get_next_state_jacobians() const
        {
            return next_state_jacobians_;
        }
        const uo_map_mx<Hessian> &uStageInternal::get_next_state_hessians() const
        {
            return next_state_hessians_;
        }
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
                from.push_back(cs::MX::vec(state));
                to.push_back(get()->get_state_syms().at(state)[k]);
            }
            for (auto &control : get()->get_controls(false))
            {
                from.push_back(cs::MX::vec(control));
                to.push_back(get()->get_control_syms().at(control)[k]);
            }
            for (auto &hybrid : get()->get_hybrids())
            {
                from.push_back(cs::MX::vec(hybrid));
                to.push_back(get()->get_hybrid_syms().at(hybrid)[k]);
            }
            for (auto &control_parameter : get()->get_control_parameters())
            {
                from.push_back(cs::MX::vec(control_parameter));
                to.push_back(get()->get_control_parameter_syms().at(control_parameter)[k]);
            }
            // for(auto& global_parameter: get()->get_global_parameters())
            // {
            //     from.push_back(global_parameter);
            //     to.push_back(global_parameter);
            // }
            return cs::MX::substitute(std::vector<cs::MX>{expr}, from, to)[0];
        }
        std::pair<std::vector<int>, cs::MX> uStage::sample(const cs::MX &expr) const
        {
            // check if expr is a column vector
            if (expr.size2() > 1)
                throw std::runtime_error("sample: expr must be a column vector");
            auto vars = cs::MX::symvar(expr);
            auto samples_vec = std::vector<cs::MX>();
            auto reti = std::vector<int>();
            if (std::all_of(vars.begin(), vars.end(), [this](const cs::MX &var)
                            { return has_variable(var); }))
            {
                for (int k = 0; k < K(); k++)
                {
                    samples_vec.push_back(eval_at_control(expr, k));
                    reti.push_back(k);
                }
            }
            return {reti, cs::MX::horzcat(samples_vec)};
        }
        const std::vector<cs::MX> &uStage::get_objective_terms() const
        {
            return get()->get_objective_terms();
        }
        std::shared_ptr<uStageInternal> uStage::get_internal() const
        {
            return static_cast<std::shared_ptr<uStageInternal>>(*this);
        }
        const std::vector<cs::MX> uStage::get_states(bool incl_auto, const std::shared_ptr<uStageInternal> &prev) const
        {
            return get()->get_states(incl_auto, prev);
        }
        const std::vector<cs::MX> uStage::get_controls(bool incl_auto, const std::shared_ptr<uStageInternal> &prev) const
        {
            return get()->get_controls(incl_auto, prev);
        }
        const std::vector<cs::MX> &uStage::get_hybrids() const
        {
            return get()->get_hybrids();
        }
        const std::vector<cs::MX> &uStage::get_control_parameters() const
        {
            return get()->get_control_parameters();
        }
        const std::vector<cs::MX> &uStage::get_global_parameters() const
        {
            return get()->get_global_parameters();
        }
        std::vector<cs::MX> uStage::all_vars() const
        {
            std::vector<cs::MX> ret;
            // get all control syms
            for (auto &control : get()->get_controls(false))
                ret.push_back(control);
            // get all state syms
            for (auto &state : get()->get_states(false))
                ret.push_back(state);
            // get all hybrid syms
            for (auto &hybrid : get()->get_hybrids())
                ret.push_back(hybrid);
            // get all control parameter syms
            for (auto &control_parameter : get()->get_control_parameters())
                ret.push_back(control_parameter);
            return ret;
        }
        void uStage::register_state(const std::vector<cs::MX> &states)
        {
            get()->register_state(states);
        };
        void uStage::register_control(const std::vector<cs::MX> &controls)
        {
            get()->register_control(controls);
        };
        void uStage::register_hybrid(const std::vector<cs::MX> &hybrids)
        {
            get()->register_hybrid(hybrids);
        };
        void uStage::register_control_parameter(const std::vector<cs::MX> &control_parameters)
        {
            uStage::get()->register_control_parameter(control_parameters);
        };
        void uStage::register_global_parameter(const std::vector<cs::MX> &global_parameters)
        {
            uStage::get()->register_global_parameter(global_parameters);
        };
        std::shared_ptr<FatropuStageEvalCasadi> uStage::get_evaluator(const std::shared_ptr<const uStageInternal> &prev, const std::shared_ptr<const uStageInternal> &next, const std::vector<cs::MX> &global_parameter_syms, const cs::Dict &opts, CasadiJitCache &cache) const
        {
            return std::make_shared<FatropuStageEvalCasadi>(uStageQuantities::create(this->get_internal(), prev, next, global_parameter_syms), opts, cache);
        }
        uStage uStage::clone() const
        {
            auto ret = uStage(*(this->get()));
            ret.get()->reset_evaluation_syms();
            ret.get()->cloned_from_ = *this;
            return ret;
        }

    }
}