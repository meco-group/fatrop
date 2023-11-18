#include "stage.hpp"
namespace fatrop
{
    namespace spectool
    {
        Stage::Stage(Ocp &ocp, const int K):std::shared_ptr<StageInternal>(std::make_shared<StageInternal>())
        {
            if (K < 2)
            {
                throw std::runtime_error("K must be greater than 1");
            }
            get()->initial = std::make_unique<uStage>(ocp.new_ustage());
            get()->middle = std::make_unique<uStage>(ocp.new_ustage(K - 2));
            get()->terminal = std::make_unique<uStage>(ocp.new_ustage());
        }

        void Stage::set_next(const cs::MX &state, const cs::MX &next_state)
        {
            at_t0().set_next(state, next_state);
            at_mid().set_next(state, next_state);
            // add the variables to the terminal stage
            at_tf().get_internal()->add_variables({state});
        }
        cs::MX Stage::sample(const cs::MX &expr)
        {
            std::vector<cs::MX> ret;
            ret.push_back(at_t0().sample(expr));
            ret.push_back(at_mid().sample(expr));
            ret.push_back(at_tf().sample(expr));
            return cs::MX::horzcat(ret);
        }
        uStage &Stage::at_t0() { return *get()->initial; }
        uStage &Stage::at_tf() { return *get()->terminal; }
        uStage &Stage::at_mid() { return *get()->middle; }
    }
}