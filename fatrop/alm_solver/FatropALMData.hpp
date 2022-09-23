#ifndef FATROPALMDATAINCLUDED
#define FATROPALMDATAINCLUDED
#include "blasfeo_wrapper/LinearAlgebraBlasfeo.hpp"
#include "templates/NLPAlg.hpp"
#include <algorithm>
using namespace std;
namespace fatrop
{
    class FatropALMData
    {
    public:
        FatropALMData(int n_ineqs) :n_ineqs_(n_ineqs), inequalitymem(n_ineqs, 5), lower_bounds(inequalitymem[0]), upper_bounds(inequalitymem[1]), ineq_curr(inequalitymem[2]), auglags_Lcurr(inequalitymem[3]), auglags_Ucurr(inequalitymem[4])
        {
        };
        int UpdateLags(double penalty_curr)
        {
            double* ineq_curr_p = ((VEC*) ineq_curr)->pa;
            double* lower_p = ((VEC*) lower_bounds)->pa;
            double* upper_p = ((VEC*) upper_bounds)->pa;
            double* auglagsL_p = ((VEC*) auglags_Lcurr)->pa;
            double* auglagsU_p = ((VEC*) auglags_Ucurr)->pa;
            double max_dist = 0.0;
            for(int i = 0; i<n_ineqs_; i ++)
            {
                double dist_lowi = ineq_curr_p[i]  - lower_p[i];
                double dist_upperi = upper_p[i] - ineq_curr_p[i];
                auglagsL_p[i] = std::max(0.0, auglagsL_p[i]- penalty_curr*dist_lowi);
                auglagsU_p[i] = std::max(0.0, auglagsU_p[i]- penalty_curr*dist_upperi);
                max_dist = std::max(std::max(max_dist, -dist_lowi), -dist_upperi);
            }
            cout << "maximum violation " << max_dist << endl;
            return 0;
        }
        const int n_ineqs_;
        FatropMemoryVecBF inequalitymem;
        FatropVecBF lower_bounds;
        FatropVecBF upper_bounds;
        FatropVecBF ineq_curr;
        FatropVecBF auglags_Lcurr;
        FatropVecBF auglags_Ucurr;
    };
}
#endif // FATROPALMDATAINCLUDED