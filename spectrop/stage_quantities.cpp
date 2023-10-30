#include "stage_quanitities.hpp"
#include "casadi_utilities.hpp"
#include "constraint_helper.hpp"
#include "stage.hpp"
#include "ocp.hpp"
namespace fatrop
{
    namespace spectrop
    {
        StageQuantities StageQuantities::create(const std::shared_ptr<StageInternal> &stage)
        {
            StageQuantities ret;
            ret.x = cs::MX::veccat(stage->get_automatics_states());
            ret.p_stage = cs::MX::veccat(stage->get_control_parameters());
            ret.u = cs::MX::veccat(stage->get_controls());
            ret.p_global = cs::MX::veccat(stage->get_ocp()->get_global_parameter_syms());
            ret.nx = ret.x.size1();
            ret.nu = ret.u.size1();
            ret.np_stage = ret.p_stage.size1();
            ret.np_global = ret.p_global.size1();
            ret.K = stage->K();
            if (stage->get_next_stage())
            {
                std::vector<cs::MX> from;
                std::vector<cs::MX> to;
                for (auto &[from_, to_] : stage->get_next_states())
                {
                    from.push_back(from_);
                    to.push_back(to_);
                }
                try
                {
                    auto x_next_vec = cs::MX::veccat(stage->get_next_stage()->get_states());
                    ret.x_next = cs::MX::veccat(cs::MX::substitute({stage->get_next_stage()->get_states()}, from, to));
                }
                catch (const std::exception &e)
                {
                    std::cerr << e.what() << '\n';
                    throw std::runtime_error("Did you set_next for every state?");
                }
                ret.nxp1 = ret.x_next.size1();
            }
            else
            {
                ret.nxp1 = 0;
            }
            ConstraintHelper::process(stage->get_constraints(),
                                      ret.lb,
                                      ret.ub,
                                      ret.g_ineq,
                                      ret.g);
            ret.ng_eq = ret.g.size1();
            ret.ng_ineq = ret.g_ineq.size1();

            ret.L = 0;
            for (auto &term : stage->get_objective_terms())
            {
                ret.L += term;
            }
            return ret;
        }
    }
}