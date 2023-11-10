#include "ustage_quanitities.hpp"
#include "casadi_utilities.hpp"
#include "constraint_helper.hpp"
#include "ustage.hpp"
#include "ocp.hpp"
namespace fatrop
{
    namespace spectrop
    {
        uStageQuantities uStageQuantities::create(const std::shared_ptr<uStageInternal> &ustage)
        {
            uStageQuantities ret;
            ret.x = cs::MX::veccat(ustage->get_states());
            ret.p_stage = cs::MX::veccat(ustage->get_control_parameters());
            ret.u = cs::MX::veccat(ustage->get_controls());
            ret.p_global = cs::MX::veccat(ustage->get_ocp()->get_global_parameter_syms());
            // ret.nx = ret.x.size1();
            // ret.nu = ret.u.size1();
            // ret.np_stage = ret.p_stage.size1();
            // ret.np_global = ret.p_global.size1();
            ret.K = ustage->K();
            if (ustage->get_next_ustage())
            {
                std::vector<cs::MX> from;
                std::vector<cs::MX> to;
                for (auto &[from_, to_] : ustage->get_next_states())
                {
                    from.push_back(from_);
                    to.push_back(to_);
                }
                // check if every state of the next stage is defined
                for(auto& state : ustage->get_next_ustage()->get_states(false))
                {
                    if(ustage->get_next_states().find(state) == ustage->get_next_states().end())
                    {
                        throw std::runtime_error("Did you set_next for every state?");
                    }
                }
                try
                {
                    auto x_next_vec = cs::MX::veccat(ustage->get_next_ustage()->get_states());
                    ret.x_next = cs::MX::veccat(cs::MX::substitute({ustage->get_next_ustage()->get_states()}, from, to));
                }
                catch (const std::exception &e)
                {
                    std::cerr << e.what() << '\n';
                    throw std::runtime_error("Did you set_next for every state?");
                }
                // ret.nxp1 = ret.x_next.size1();
            }
            else
            {
                // ret.nxp1 = 0;
                ret.x_next = cs::MX::zeros(0, 1);
            }
            ConstraintHelper::process(ustage->get_constraints(),
                                      ret.lb,
                                      ret.ub,
                                      ret.g_ineq,
                                      ret.g);
            // ret.ng_eq = ret.g.size1();
            // ret.ng_ineq = ret.g_ineq.size1();

            ret.L = 0;
            for (auto &term : ustage->get_objective_terms())
            {
                ret.L += term;
            }
            return ret;
        }
    }
}