#pragma once
#include "ocp.hpp"
#include "ustage.hpp"
#include <memory>

namespace fatrop
{
    namespace spectool
    {
        enum class at
        {
            t0,
            mid,
            tf
        };
        class Stage; // forward declaration
        class StageInternal
        {
            friend class Stage;

            std::unique_ptr<uStage> initial;
            std::unique_ptr<uStage> middle;
            std::unique_ptr<uStage> terminal;
        };
        class Stage : public std::shared_ptr<StageInternal>
        {
            friend class Ocp;

        private:
            Stage(Ocp &ocp, const int K);

        public:
            void set_next(const cs::MX &state, const cs::MX &next_state, const Jacobian & jacobian = Jacobian{}, const Hessian & hessian = Hessian{});
            // void set_next(const cs::MX &state, const cs::MX &next_state);
            std::pair<std::vector<int>, cs::MX> sample(const cs::MX &expr) const;
            template <class... args>
            void subject_to(const cs::MX &expr, args... argss)
            {
                apply_at(static_cast<void (uStage::*)(const cs::MX&)>(&uStage::subject_to), expr, argss...);
            }
            template <class... args>
            void add_objective(const cs::MX &expr, args... argss)
            {
                apply_at(&uStage::add_objective, expr, argss...);
            }
            uStage &at_t0() const;
            uStage &at_tf() const;
            uStage &at_mid() const;
            cs::MX at_t0(const cs::MX &expr) const { return at_t0().at_t0(expr); };
            cs::MX at_tf(const cs::MX &expr) const { return at_tf().at_t0(expr); };

            template <class F>
            void apply_at_single(F f, const cs::MX &expr, const at &type)
            {
                switch (type)
                {
                case at::t0:
                    (at_t0().*f)(expr);
                    break;
                case at::mid:
                    (at_mid().*f)(expr);
                    break;
                case at::tf:
                    (at_tf().*f)(expr);
                    break;
                default:
                    break;
                }
            }
            template <typename F, class... args>
            void apply_at(F f, const cs::MX& expr, args... argss)
            {
                // count the number of arguments
                int n = sizeof...(args);
                // create array of arguments
                at arr[] = {argss...};
                // apply for every argument
                for (int i = 0; i < n; i++)
                {
                    apply_at_single(f, expr, arr[i]);
                }
            }

        private:
        };

    }

}