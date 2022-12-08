#ifndef FATROPALMALGINCLUDED
#define FATROPALMALGINCLUDED
using namespace std;
#include "solver/FatropAlg.hpp"
#include "templates/FatropNLPAL.hpp"
#include "FatropALMData.hpp"
#include <memory>
#include <templates/FatropApplication.hpp>
using namespace std;
namespace fatrop
{
    class FatropALMAlg : public FatropApplication
    {
    public:
        FatropALMAlg(
            const shared_ptr<FatropNLPAL> &fatropnlpal,
            const shared_ptr<FatropData> &fatropdata,
            const shared_ptr<FatropParams> &fatropparams,
            const shared_ptr<Filter> &filter,
            const shared_ptr<LineSearch> &linesearch,
            const shared_ptr<Journaller> &journaller) : fatropnlpal_(fatropnlpal), innersolver_(fatropnlpal,
                                                                                                fatropdata,
                                                                                                fatropparams,
                                                                                                filter,
                                                                                                linesearch,
                                                                                                journaller),
                                                        almdata_(fatropnlpal->GetNOIneqs(), fatropnlpal->GetNLPDims().nvars)
        {
            fatropdata->x_initial.copyto(almdata_.initial_x);
        }
        void Initialize() override{};
        void Reset() override{};
        void SetBounds(const vector<double> &lower_boundsin, const vector<double> &upper_boundsin) override;
        void SetInitial(const vector<double> &initial) override;
        void GetSolution(vector<double> &sol) override;
        int Optimize() override;
        void WarmStart() override
        {
            innersolver_.fatropdata_->x_curr.copyto(almdata_.initial_x);
        }
        const shared_ptr<FatropNLPAL> fatropnlpal_;
        FatropAlg innersolver_;
        FatropALMData almdata_;
    };
} // namespace fatrop
#endif // FATROPALMALGINCLUDED