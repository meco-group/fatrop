#include "ustage_quantities.hpp"
#include "fatrop/spectool/auxiliary/casadi_utilities.hpp"
#include "fatrop/spectool/auxiliary/constraint_helper.hpp"
#include "ustage.hpp"
#include "ocp.hpp"
namespace fatrop
{
    namespace spectool
    {
        uStageQuantities uStageQuantities::create(const std::shared_ptr<const uStageInternal> &ustage, const std::shared_ptr<const uStageInternal> &prev, const std::shared_ptr<const uStageInternal> &next)
        {
            uStageQuantities ret;
            ret.x = cs::MX::veccat(ustage->get_states(true, prev));
            ret.p_stage = cs::MX::veccat(ustage->get_control_parameters());
            ret.u = cs::MX::veccat(ustage->get_controls(true, prev));
            ret.p_global = cs::MX::veccat(ustage->get_ocp()->get_global_parameter_syms());
            ret.K = ustage->K();
            // std::cout << "number of hybrids that are states " << ustage->get_hybrids_states(prev).size() << std::endl;
            // std::cout << "number of hybrids that are controls " << ustage->get_hybrids_controls(prev).size() << std::endl;
            if (next)
            {
                std::vector<cs::MX> from;
                std::vector<cs::MX> to;
                for (auto &[from_, to_] : ustage->get_next_states())
                {
                    from.push_back(from_);
                    to.push_back(to_);
                }
                // check if every state of the next stage is defined
                for(auto& state : next->get_states(false))
                {
                    if(ustage->get_next_states().find(state) == ustage->get_next_states().end())
                    {
                        throw std::runtime_error("Did you set_next for every state?");
                    }
                }
                try
                {
                    auto x_next_vec = cs::MX::veccat(next->get_states(true, ustage));
                    ret.x_next = cs::MX::veccat(cs::MX::substitute({next->get_states(true, ustage)}, from, to));
                }
                catch (const std::exception &e)
                {
                    std::cerr << e.what() << '\n';
                    throw std::runtime_error("Did you set_next for every state?");
                }
                if(ustage->K()>1)
                {
                    // check if x_next is same as x
                    auto x_next_vec = next->get_states(true, ustage);
                    auto x_vec = ustage->get_states(true, prev);
                    // check if size is same
                    if(x_next_vec.size() != x_vec.size())
                    {
                        throw std::runtime_error("x_next and x must have the same size");
                    }
                    // iterate over elements and check if the same
                    for(int i = 0; i < x_next_vec.size(); i++)
                    {
                        if (x_next_vec[i].get() != x_vec[i].get())
                        {
                            throw std::runtime_error("x_next and x must be the same");
                        }
                    }
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