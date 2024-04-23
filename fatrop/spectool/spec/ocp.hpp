#pragma once
#include <casadi/casadi.hpp>
#include <string>
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
#include "ustage.hpp"
#include "fatrop/spectool/solver_interfaces/solver.hpp"
#include "stage.hpp"
#include <utility>
// #include "ustage_eval_casadi.hpp"

namespace fatrop
{
    namespace spectool
    {
        class Stage;
        namespace cs = casadi;
        class Ocp;
        // class Stage;
        class OcpInternal
        {
        public:
            const uo_set_mx &get_states();
            const uo_set_mx &get_controls();
            const uo_set_mx &get_global_parameters();
            const uo_set_mx &get_control_parameters();
            const std::vector<cs::MX> &get_global_parameter_syms()
            {
                return global_parammeter_syms_;
            }

        protected:
            friend class Ocp;
            friend class uStageInternal;
            uo_set_mx states_;
            uo_set_mx controls_;
            uo_set_mx hybrids_;
            uo_set_mx global_parameters_;
            uo_set_mx control_parameters_;
            std::vector<cs::MX> global_parammeter_syms_;
            bool is_state(const cs::MX &var);
            bool is_control(const cs::MX &var);
            bool is_hybrid(const cs::MX &var);
            bool is_global_parameter(const cs::MX &var);
            bool is_control_parameter(const cs::MX &var);
            std::vector<std::pair<cs::MX, cs::MX>> initial_values;
            std::vector<std::pair<cs::MX, cs::MX>> parameter_values;
            std::string solver_name = "";
            std::vector<cs::MX> multi_ustage_constraints;
            std::shared_ptr<SolverInterface> solver_ptr;

        private:
            void add_to_ordering(const cs::MX &var);
            std::vector<cs::MX> order_vars(const std::vector<cs::MX> &vars);
            uo_map_mx<size_t> ordering_;
        };
        class Ocp : private std::shared_ptr<OcpInternal>
        {
        public:
            Ocp() : std::shared_ptr<OcpInternal>(new OcpInternal())
            {
            }
            cs::MX state(const int m = 1, const int n = 1);
            cs::MX control(const int m = 1, const int n = 1);
            cs::MX hybrid(const int m = 1, const int n = 1);
            cs::MX parameter(const int m = 1, const int n = 1, const std::string &grid = "global");
            std::pair<std::vector<int>,cs::MX> sample(const cs::MX &expr) const;
            uStage new_ustage(const int K = 1);
            Stage new_stage(const int K);
            void add_ustage(const uStage &ustage);
            const std::vector<uStage> &get_ustages() const;
            const std::vector<cs::MX> &get_global_parameters() const
            {
                return get()->global_parammeter_syms_;
            }
            const std::vector<cs::MX> &get_multi_ustage_constraints() const
            {
                return get()->multi_ustage_constraints;
            }
            cs::Function to_function(const std::string& name, const std::vector<cs::MX> &in, const std::vector<cs::MX> &out) const;
            cs::MX at_t0(const cs::MX &expr) const { return ustages_.front().at_t0(expr); };
            cs::MX at_tf(const cs::MX &expr) const { return ustages_.back().at_tf(expr); };
            uStage at_t0() const { return ustages_.front(); };
            uStage at_tf() const { return ustages_.back(); };
            void set_initial(const cs::MX &var, const cs::MX &value);
            std::vector<cs::MX> eval_at_initial(const std::vector<cs::MX> &expr) const;
            cs::MX all_variables() const;
            void solver(const std::string &name, const cs::Dict &function_opts = casadi::Dict(), const cs::Dict &solver_opts = casadi::Dict()){get()->solver_name = name; solver_opts_ = solver_opts; function_opts_ = function_opts;};
            std::shared_ptr<SolverInterface> get_solver() const {return get()->solver_ptr;};
            void subject_to(const cs::MX &expr);

        protected:
            std::vector<uStage> ustages_;
            cs::Dict function_opts_;
            cs::Dict solver_opts_;

        };
    } // namespace spectrop
} // namespace fatrop