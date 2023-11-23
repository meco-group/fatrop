#include "stage.hpp"
namespace fatrop
{
    namespace spectool
    {
        Stage::Stage(Ocp &ocp, const int K) : std::shared_ptr<StageInternal>(std::make_shared<StageInternal>())
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
        std::pair<cs::MX, std::vector<int>> Stage::sample(const cs::MX &expr) const
        {
            if (expr.size2() != 1)
                return {cs::MX(), {}}; // return empty matrix if input is not a column vector
            auto ret = std::vector<cs::MX>();
            auto vars = cs::MX::symvar(expr);
            auto reti = std::vector<int>();
            int k = 0;
            auto ustages_ = std::array<const uStage*, 3>{&at_t0(), &at_mid(), &at_tf()};
            for (auto i = 0; i < ustages_.size(); ++i)
            {
                const auto &ustage = ustages_[i];
                {
                    auto sample_ = (*ustage).sample(expr);
                    if(sample_.first.size2() == 0)
                        continue;
                    // add samples to ret
                    ret.push_back(sample_.first);
                    // add indices to reti
                    for (const auto &idx : sample_.second)
                        reti.push_back(k + idx);
                }
                k += (*ustage).K();
            }
            return {cs::MX::horzcat(ret), reti};
        }
        uStage &Stage::at_t0() const { return *get()->initial; }
        uStage &Stage::at_tf() const { return *get()->terminal; }
        uStage &Stage::at_mid() const { return *get()->middle; }
    }
}