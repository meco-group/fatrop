// Copyright (c) 2026 Lander Vanroye, KU Leuven
#ifndef FATROP_IP_MULTIPLIER_RECALCULATOR_HPP
#define FATROP_IP_MULTIPLIER_RECALCULATOR_HPP

#include "fatrop/common/options.hpp"
#include "fatrop/ip_algorithm/ip_data.hpp"
#include "fatrop/ip_algorithm/ip_eq_mult_initializer.hpp"
#include <cmath>
#include <memory>
#include <stdexcept>

namespace fatrop
{
    // Recalculate equality multipliers after accepting a sufficiently feasible point.
    // The default follows Ipopt: enabled for limited-memory Hessians only.
    template <typename ProblemType> class IpMultiplierRecalculator
    {
    public:
        IpMultiplierRecalculator(std::shared_ptr<IpData<ProblemType>> data,
                                std::shared_ptr<IpEqMultInitializerBase> calculator)
            : data_(std::move(data)), calculator_(std::move(calculator)) {}

        void register_options(OptionRegistry &registry)
        {
            registry.register_option("recalc_y", &IpMultiplierRecalculator::set_enabled, this);
            registry.register_option("recalc_y_feas_tol", &IpMultiplierRecalculator::set_feas_tol,
                                     this);
        }

        void set_enabled(const bool &enabled)
        {
            enabled_ = enabled;
            explicitly_set_ = true;
        }

        void set_feas_tol(const Scalar &tol)
        {
            if (!(tol > 0.) || !std::isfinite(tol))
                throw std::invalid_argument("recalc_y_feas_tol must be positive and finite");
            feas_tol_ = tol;
        }

        void after_accepted_step()
        {
            if (!(explicitly_set_ ? enabled_ : data_->hessian_evaluator().uses_history())) return;
            auto &iterate = data_->current_iterate();
            if (iterate.info().number_of_eq_constraints == 0) return;
            if (iterate.constr_viol_inf() <= feas_tol_) calculator_->recalculate_eq_mult();
        }

    private:
        std::shared_ptr<IpData<ProblemType>> data_;
        std::shared_ptr<IpEqMultInitializerBase> calculator_;
        Scalar feas_tol_ = 1e-6;
        bool enabled_ = false;
        bool explicitly_set_ = false;
    };
}
#endif
