//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_timings_hpp__
#define __fatrop_ip_algorithm_ip_timings_hpp__

#include "fatrop/common/printing.hpp"
#include "fatrop/common/timing.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace fatrop
{
    struct IpTimingStatistics
    {
        // High-level timers
        Timer full_algorithm;
        Timer initialization;
        Timer compute_search_dir;
        // NLP-eval timers (driven from ip_nlp_orig.hxx)
        Timer eval_objective;
        Timer eval_gradient;
        Timer eval_constraint_violation;
        Timer eval_hessian;
        Timer eval_jacobian;
        // Iterate-level computation timers (driven from ip_iterate.hxx).
        // Each timer measures the exclusive cost of its method: time spent
        // inside any nested timed scope (another iterate timer or an
        // eval_* timer) is excluded via the active_stack_ mechanism below.
        Timer compute_obj_value;
        Timer compute_obj_gradient;
        Timer compute_constr_viol;
        Timer compute_hessian;
        Timer compute_jacobian;
        Timer compute_barrier_value;
        Timer compute_barrier_gradient;
        Timer compute_delta_lower;
        Timer compute_delta_upper;
        Timer compute_dual_infeasibility_x;
        Timer compute_dual_infeasibility_s;
        Timer compute_complementarity_l;
        Timer compute_complementarity_u;
        Timer compute_relaxed_compl_l;
        Timer compute_relaxed_compl_u;
        Timer compute_primal_damping;
        Timer compute_linear_decrease_obj;
        Timer compute_linear_decrease_bar;
        Timer compute_e_mu;
        Timer compute_max_step_primal;
        Timer compute_max_step_dual;
        Timer compute_modify_dual_bounds;

        // Currently active scoped timer (or nullptr if none).  When a
        // ScopedTimer is constructed it snapshots this pointer into its
        // own field, pauses it (if non-null), and overwrites this field
        // with itself.  On destruction the snapshot is restored.  This
        // gives nested-timer exclusion without any heap allocation.
        Timer *active_ = nullptr;

        inline double compute_rest_time() const;
        inline double compute_function_evaluation() const;
        inline double compute_iterate_quantities() const;
        inline double compute_fatrop() const;
        inline void reset();
        inline void pause();
        friend inline std::ostream &operator<<(std::ostream &os, const IpTimingStatistics &timings);
    };

    // RAII helper that pauses the currently active timer (if any),
    // starts a new one, and on destruction resumes the previously active
    // timer.  The "previous" pointer is stored in this ScopedTimer's own
    // stack frame, so the chain of nested ScopedTimers behaves like a
    // linked list threaded through the call stack -- no heap allocation.
    // Both NLP-eval timers (in ip_nlp_orig.hxx) and iterate-level timers
    // (in ip_iterate.hxx) use this so they participate in the same
    // exclusion chain.
    struct ScopedTimer
    {
        ScopedTimer(Timer &t, IpTimingStatistics &s)
            : timer_(&t), prev_(s.active_), stats_(&s)
        {
            if (prev_)
                prev_->pause();
            timer_->start();
            stats_->active_ = timer_;
        }
        ~ScopedTimer()
        {
            timer_->pause();
            stats_->active_ = prev_;
            if (prev_)
                prev_->start();
        }
        ScopedTimer(const ScopedTimer &) = delete;
        ScopedTimer &operator=(const ScopedTimer &) = delete;

      private:
        Timer *timer_;
        Timer *prev_;
        IpTimingStatistics *stats_;
    };

    void IpTimingStatistics::reset()
    {
        full_algorithm.reset();
        initialization.reset();
        compute_search_dir.reset();
        eval_objective.reset();
        eval_gradient.reset();
        eval_constraint_violation.reset();
        eval_hessian.reset();
        eval_jacobian.reset();
        compute_obj_value.reset();
        compute_obj_gradient.reset();
        compute_constr_viol.reset();
        compute_hessian.reset();
        compute_jacobian.reset();
        compute_barrier_value.reset();
        compute_barrier_gradient.reset();
        compute_delta_lower.reset();
        compute_delta_upper.reset();
        compute_dual_infeasibility_x.reset();
        compute_dual_infeasibility_s.reset();
        compute_complementarity_l.reset();
        compute_complementarity_u.reset();
        compute_relaxed_compl_l.reset();
        compute_relaxed_compl_u.reset();
        compute_primal_damping.reset();
        compute_linear_decrease_obj.reset();
        compute_linear_decrease_bar.reset();
        compute_e_mu.reset();
        compute_max_step_primal.reset();
        compute_max_step_dual.reset();
        compute_modify_dual_bounds.reset();
        active_ = nullptr;
    }

    void IpTimingStatistics::pause()
    {
        full_algorithm.pause();
        initialization.pause();
        compute_search_dir.pause();
        eval_objective.pause();
        eval_gradient.pause();
        eval_constraint_violation.pause();
        eval_hessian.pause();
        eval_jacobian.pause();
        compute_obj_value.pause();
        compute_obj_gradient.pause();
        compute_constr_viol.pause();
        compute_hessian.pause();
        compute_jacobian.pause();
        compute_barrier_value.pause();
        compute_barrier_gradient.pause();
        compute_delta_lower.pause();
        compute_delta_upper.pause();
        compute_dual_infeasibility_x.pause();
        compute_dual_infeasibility_s.pause();
        compute_complementarity_l.pause();
        compute_complementarity_u.pause();
        compute_relaxed_compl_l.pause();
        compute_relaxed_compl_u.pause();
        compute_primal_damping.pause();
        compute_linear_decrease_obj.pause();
        compute_linear_decrease_bar.pause();
        compute_e_mu.pause();
        compute_max_step_primal.pause();
        compute_max_step_dual.pause();
        compute_modify_dual_bounds.pause();
    }

    double IpTimingStatistics::compute_rest_time() const
    {
        double ret = full_algorithm.elapsed();
        ret -= initialization.elapsed();
        ret -= compute_search_dir.elapsed();
        ret -= eval_objective.elapsed();
        ret -= eval_gradient.elapsed();
        ret -= eval_constraint_violation.elapsed();
        ret -= eval_hessian.elapsed();
        ret -= eval_jacobian.elapsed();
        ret -= compute_iterate_quantities();
        return ret;
    }

    double IpTimingStatistics::compute_function_evaluation() const
    {
        double ret = 0.;
        ret += eval_objective.elapsed();
        ret += eval_gradient.elapsed();
        ret += eval_constraint_violation.elapsed();
        ret += eval_hessian.elapsed();
        ret += eval_jacobian.elapsed();
        return ret;
    }

    double IpTimingStatistics::compute_iterate_quantities() const
    {
        double ret = 0.;
        ret += compute_obj_value.elapsed();
        ret += compute_obj_gradient.elapsed();
        ret += compute_constr_viol.elapsed();
        ret += compute_hessian.elapsed();
        ret += compute_jacobian.elapsed();
        ret += compute_barrier_value.elapsed();
        ret += compute_barrier_gradient.elapsed();
        ret += compute_delta_lower.elapsed();
        ret += compute_delta_upper.elapsed();
        ret += compute_dual_infeasibility_x.elapsed();
        ret += compute_dual_infeasibility_s.elapsed();
        ret += compute_complementarity_l.elapsed();
        ret += compute_complementarity_u.elapsed();
        ret += compute_relaxed_compl_l.elapsed();
        ret += compute_relaxed_compl_u.elapsed();
        ret += compute_primal_damping.elapsed();
        ret += compute_linear_decrease_obj.elapsed();
        ret += compute_linear_decrease_bar.elapsed();
        ret += compute_e_mu.elapsed();
        ret += compute_max_step_primal.elapsed();
        ret += compute_max_step_dual.elapsed();
        ret += compute_modify_dual_bounds.elapsed();
        return ret;
    }

    double IpTimingStatistics::compute_fatrop() const
    {
        return full_algorithm.elapsed() - compute_function_evaluation();
    }

    std::ostream &operator<<(std::ostream &os, const IpTimingStatistics &timings)
    {
        const int name_w = 36;
        const int time_w = 12;
        const int pct_w = 8;
        const int line_w = name_w + time_w + pct_w;
        const double total = timings.full_algorithm.elapsed();
        auto fmt_time = [](double v) {
            std::ostringstream s;
            s << std::fixed << std::setprecision(6) << v;
            return s.str();
        };
        auto fmt_pct = [&](double v) {
            std::ostringstream s;
            if (total > 0.)
                s << std::fixed << std::setprecision(1) << (100. * v / total) << "%";
            else
                s << "-";
            return s.str();
        };
        auto row = [&os, name_w, time_w](const char *name, double v) {
            os << std::left << std::setw(name_w) << name << std::right << std::setw(time_w)
               << std::fixed << std::setprecision(6) << v << std::endl;
        };
        auto row_indent = [&os, name_w, time_w](const char *name, double v) {
            std::string indented = std::string("  ") + name;
            os << std::left << std::setw(name_w) << indented << std::right << std::setw(time_w)
               << std::fixed << std::setprecision(6) << v << std::endl;
        };
        auto row_pct = [&](const char *name, double v) {
            os << std::left << std::setw(name_w) << name << std::right << std::setw(time_w)
               << fmt_time(v) << std::right << std::setw(pct_w) << fmt_pct(v) << std::endl;
        };
        os << std::left << std::setw(name_w) << "operation" << std::right << std::setw(time_w)
           << "time (s)" << std::right << std::setw(pct_w) << "share" << std::endl;
        os << std::string(line_w, '-') << std::endl;
        row_pct("fatrop", timings.compute_fatrop());
        row_indent("initialization", timings.initialization.elapsed());
        row_indent("compute search dir", timings.compute_search_dir.elapsed());
        row_indent("compute obj value", timings.compute_obj_value.elapsed());
        row_indent("compute obj gradient", timings.compute_obj_gradient.elapsed());
        row_indent("compute constr viol", timings.compute_constr_viol.elapsed());
        row_indent("compute hessian", timings.compute_hessian.elapsed());
        row_indent("compute jacobian", timings.compute_jacobian.elapsed());
        row_indent("compute barrier value", timings.compute_barrier_value.elapsed());
        row_indent("compute barrier gradient", timings.compute_barrier_gradient.elapsed());
        row_indent("compute delta lower", timings.compute_delta_lower.elapsed());
        row_indent("compute delta upper", timings.compute_delta_upper.elapsed());
        row_indent("compute dual infeasibility x", timings.compute_dual_infeasibility_x.elapsed());
        row_indent("compute dual infeasibility s", timings.compute_dual_infeasibility_s.elapsed());
        row_indent("compute complementarity l", timings.compute_complementarity_l.elapsed());
        row_indent("compute complementarity u", timings.compute_complementarity_u.elapsed());
        row_indent("compute relaxed compl l", timings.compute_relaxed_compl_l.elapsed());
        row_indent("compute relaxed compl u", timings.compute_relaxed_compl_u.elapsed());
        row_indent("compute primal damping", timings.compute_primal_damping.elapsed());
        row_indent("compute linear decrease obj", timings.compute_linear_decrease_obj.elapsed());
        row_indent("compute linear decrease bar", timings.compute_linear_decrease_bar.elapsed());
        row_indent("compute e_mu", timings.compute_e_mu.elapsed());
        row_indent("compute max step primal", timings.compute_max_step_primal.elapsed());
        row_indent("compute max step dual", timings.compute_max_step_dual.elapsed());
        row_indent("compute modify dual bounds", timings.compute_modify_dual_bounds.elapsed());
        row_indent("rest", timings.compute_rest_time());
        row_pct("function evaluation", timings.compute_function_evaluation());
        row_indent("eval objective", timings.eval_objective.elapsed());
        row_indent("eval gradient", timings.eval_gradient.elapsed());
        row_indent("eval constraint violation", timings.eval_constraint_violation.elapsed());
        row_indent("eval hessian", timings.eval_hessian.elapsed());
        row_indent("eval jacobian", timings.eval_jacobian.elapsed());
        os << std::string(line_w, '-') << std::endl;
        row("total", total);
        return os;
    }

} // namespace fatrop

#endif // __fatrop_ip_algorithnm_ip_timings_hpp__
