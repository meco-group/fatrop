
#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include <memory>
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
#include "fatrop/spectool/function_evaluation/casadi_jit_cache.hpp"
#include "ustage_eval_abstract.hpp"
#include "ustage_eval_casadi.hpp"
namespace fatrop
{
    namespace spectool
    {
        namespace cs = casadi;
        class uStage;
        class Stage;
        class OcpInternal;
        class uStageInternal
        {
        public:
            uStageInternal(int K, const std::shared_ptr<OcpInternal> &ocp) : K_(K), ocp_(ocp), auto_mode(true)
            {
            }
            uStageInternal(int K, const std::vector<cs::MX> &states, const std::vector<cs::MX> &controls, const std::vector<cs::MX> &control_parameters, const std::vector<cs::MX> &global_parameters) : K_(K), auto_mode(false)
            {
                register_state(states);
                register_control(controls);
                register_control_parameter(control_parameters);
                register_global_parameter(global_parameters);
            }
            uStageInternal(int K) : K_(K), auto_mode(false) // auto mode only possible with reference to ocp
                                    {};
            uStageInternal() : uStageInternal{1, {}, {}, {}, {}} {}; // empty ustage
            friend class uStage;
            friend class Stage;
            friend class Ocp;

        protected:
            void add_variables(const cs::MX &expr);

        public:
            const std::vector<cs::MX> &get_objective_terms() const;
            const std::vector<cs::MX> &get_constraints() const;
            const uo_map_mx<cs::MX> &get_next_states() const;
            const std::vector<cs::MX> get_states(bool include_hybrids = true, const std::shared_ptr<const uStageInternal> &prev = nullptr) const;
            const uo_map_mx<std::vector<cs::MX>> &get_state_syms() const;
            const std::vector<cs::MX> get_controls(bool include_hybrids = true, const std::shared_ptr<const uStageInternal> &prev = nullptr) const;
            const uo_map_mx<std::vector<cs::MX>> &get_control_syms() const;
            const std::vector<cs::MX> &get_hybrids() const;
            const std::vector<cs::MX> get_hybrids_states(const std::shared_ptr<const uStageInternal> &prev) const;
            const std::vector<cs::MX> get_hybrids_controls(const std::shared_ptr<const uStageInternal> &prev) const;
            const uo_map_mx<std::vector<cs::MX>> &get_hybrid_syms() const;
            const std::vector<cs::MX> &get_control_parameters() const;
            const uo_map_mx<std::vector<cs::MX>> &get_control_parameter_syms() const;
            const std::vector<cs::MX> &get_global_parameters() const;
            const std::shared_ptr<OcpInternal> get_ocp() const;
            int K() const { return K_; };

        private:
            bool has_variable(const cs::MX &var) const;
            void register_state(const cs::MX &state);
            void register_control(const cs::MX &control);
            void register_hybrid(const cs::MX &hybrid);
            void register_control_parameter(const cs::MX &control_parameter);
            void register_global_parameter(const cs::MX &global_parameter);
            void register_state(const std::vector<cs::MX> &states);
            void register_control(const std::vector<cs::MX> &controls);
            void register_hybrid(const std::vector<cs::MX> &hybrids);
            void register_control_parameter(const std::vector<cs::MX> &control_parameters);
            void register_global_parameter(const std::vector<cs::MX> &global_parameters);
            void get_hybrids(std::vector<cs::MX> &states, std::vector<cs::MX> &controls, const std::shared_ptr<const uStageInternal> &prev) const;
            void reset_evaluation_syms()
            {
                state_syms_.clear();
                control_syms_.clear();
                hybrid_syms_.clear();
                control_parameter_syms_.clear();
                register_state(states_);
                register_control(controls_);
                register_hybrid(hybrids_);
                register_control_parameter(control_parameters_);
            }
            const int K_;
            std::vector<cs::MX> objective_terms_;
            std::vector<cs::MX> constraints_;
            uo_map_mx<cs::MX> next_states_;
            std::vector<cs::MX> states_;
            uo_set_mx states_set_;
            uo_map_mx<std::vector<cs::MX>> state_syms_;
            std::vector<cs::MX> controls_;
            uo_set_mx controls_set_;
            std::vector<cs::MX> hybrids_;
            uo_set_mx hybrids_set_;
            uo_map_mx<std::vector<cs::MX>> control_syms_;
            uo_map_mx<std::vector<cs::MX>> hybrid_syms_;
            std::vector<cs::MX> control_parameters_;
            uo_set_mx control_parameters_set_;
            uo_map_mx<std::vector<cs::MX>> control_parameter_syms_;
            std::weak_ptr<OcpInternal> ocp_; // is only used for determining symbol type (control/state/hybrid/parameter) in automatic mode
            std::vector<cs::MX> global_parameters_;
            uo_set_mx global_parameters_set_;
            bool auto_mode = false;
            bool has_hybrids = false;
        };
        class uStage : private std::shared_ptr<uStageInternal>
        {
        public:
            friend class Ocp;
            template <class... Args>
            uStage(Args &&...args) : std::shared_ptr<uStageInternal>(new uStageInternal(std::forward<Args>(args)...))
            {
            }
            void subject_to(const cs::MX &constraint);
            void add_objective(const cs::MX &objective);
            void set_next(const cs::MX &state, const cs::MX &next_state);
            void set_next(const uo_map_mx<cs::MX> &next_states);
            bool has_variable(const cs::MX &var) const;
            cs::MX at_t0(const cs::MX &expr) const { return eval_at_control(expr, 0); };
            cs::MX at_tf(const cs::MX &expr) const { return eval_at_control(expr, K() - 1); };
            int K() const;
            const uo_map_mx<cs::MX> &dynamics();
            cs::MX eval_at_control(const cs::MX &expr, const int k) const;
            std::pair<std::vector<int>, cs::MX> sample(const cs::MX &expr) const;
            std::shared_ptr<uStageInternal> get_internal() const;
            const std::vector<cs::MX> &get_objective_terms() const;
            const std::vector<cs::MX> get_states(bool include_hybrids = true, const std::shared_ptr<uStageInternal> &prev = nullptr) const;
            const std::vector<cs::MX> get_controls(bool include_hybrids = true, const std::shared_ptr<uStageInternal> &prev = nullptr) const;
            const std::vector<cs::MX> &get_hybrids() const;
            const std::vector<cs::MX> &get_control_parameters() const;
            const std::vector<cs::MX> &get_global_parameters() const;
            uStage clone() const;
            void register_state(const std::vector<cs::MX> &states);
            void register_control(const std::vector<cs::MX> &controls);
            void register_hybrid(const std::vector<cs::MX> &hybrids);
            void register_control_parameter(const std::vector<cs::MX> &control_parameters);
            void register_global_parameter(const std::vector<cs::MX> &global_parameters);
            virtual std::shared_ptr<FatropuStageEvalAbstract> get_evaluator(const std::shared_ptr<const uStageInternal> &prev, const std::shared_ptr<const uStageInternal> &next, const std::vector<cs::MX> &global_parameter_syms, const cs::Dict &opts, CasadiJitCache &cache) const;
            bool operator==(const uStage &other) const { return get_internal() == other.get_internal(); };

            // private:
            //     uStage(const std::shared_ptr<uStageInternal> &internal) : std::shared_ptr<uStageInternal>(internal){};
        };
    } // namespace spectrop
} // namespace fatrop