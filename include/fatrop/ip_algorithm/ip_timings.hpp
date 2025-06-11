//
// Copyright (c) 2024 Lander Vanroye, KU Leuven
//

#ifndef __fatrop_ip_algorithm_ip_timings_hpp__
#define __fatrop_ip_algorithm_ip_timings_hpp__

#include "fatrop/common/printing.hpp"
#include "fatrop/common/timing.hpp"
#include <iomanip>
#include <iostream>

namespace fatrop
{
    struct IpTimingStatistics
    {
        Timer full_algorithm;
        Timer initialization;
        Timer compute_search_dir;
        Timer eval_objective;
        Timer eval_gradient;
        Timer eval_constraint_violation;
        Timer eval_hessian;
        Timer eval_jacobian;
        inline double compute_rest_time() const;
        inline double compute_function_evaluation() const;
        inline double compute_fatrop() const;
        inline void reset();
        inline void pause();
        friend inline std::ostream &operator<<(std::ostream &os, const IpTimingStatistics &timings);
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

    double IpTimingStatistics::compute_fatrop() const
    {
        return full_algorithm.elapsed() - compute_function_evaluation();
    }

    std::ostream &operator<<(std::ostream &os, const IpTimingStatistics &timings)
    {
        os << std::left << std::setw(30) << "operation" << std::right << std::setw(10) << "time (s)"
           << std::endl;
        os << std::string(40, '-') << std::endl;
        os << std::left << std::setw(30) << "initialization" << std::right << std::setw(10)
           << std::fixed << std::setprecision(6) << timings.initialization.elapsed() << std::endl;
        os << std::left << std::setw(30) << "compute search dir" << std::right << std::setw(10)
           << std::fixed << std::setprecision(6) << timings.compute_search_dir.elapsed() << std::endl;
        os << std::left << std::setw(30) << "eval objective" << std::right << std::setw(10)
           << std::fixed << std::setprecision(6) << timings.eval_objective.elapsed() << std::endl;
        os << std::left << std::setw(30) << "eval gradient" << std::right << std::setw(10)
           << std::fixed << std::setprecision(6) << timings.eval_gradient.elapsed() << std::endl;
        os << std::left << std::setw(30) << "eval constraint violation" << std::right
           << std::setw(10) << std::fixed << std::setprecision(6)
           << timings.eval_constraint_violation.elapsed() << std::endl;
        os << std::left << std::setw(30) << "eval hessian" << std::right << std::setw(10)
           << std::fixed << std::setprecision(6) << timings.eval_hessian.elapsed() << std::endl;
        os << std::left << std::setw(30) << "eval jacobian" << std::right << std::setw(10)
           << std::fixed << std::setprecision(6) << timings.eval_jacobian.elapsed() << std::endl;
        os << std::left << std::setw(30) << "rest" << std::right << std::setw(10) << std::fixed
           << std::setprecision(6) << timings.compute_rest_time() << std::endl;
        os << std::string(40, '-') << std::endl;
        os << std::left << std::setw(30) << "time function evaluation" << std::right
           << std::setw(10) << std::fixed << std::setprecision(6)
           << timings.compute_function_evaluation() << std::endl;
        os << std::left << std::setw(30) << "time fatrop" << std::right << std::setw(10)
           << std::fixed << std::setprecision(6) << timings.compute_fatrop() << std::endl;
        os << std::string(40, '-') << std::endl;
        os << std::left << std::setw(30) << "total " << std::right << std::setw(10)
           << std::fixed << std::setprecision(6) << timings.full_algorithm.elapsed() << std::endl;
        return os;
    }

} // namespace fatrop

#endif // __fatrop_ip_algorithnm_ip_timings_hpp__
